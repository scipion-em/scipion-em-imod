# *****************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# *****************************************************************************

import os

from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
import tomo.objects as tomoObj
from pwem.emlib.image import ImageHandler

from .. import Plugin, utils
from .protocol_base import ProtImodBase, EXT_MRCS_TS_EVEN_NAME, EXT_MRCS_TS_ODD_NAME

SCIPION_IMPORT = 0
FIXED_DOSE = 1


class ProtImodDoseFilter(ProtImodBase):
    """
    Tilt-series dose filtering based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/mtffilter.html
    """

    _label = 'Dose filter'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('initialDose',
                      params.FloatParam,
                      default=0.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Initial dose (e/sq. Å)',
                      help='Dose applied before any of the images in the '
                           'input file were taken; this value will be '
                           'added to all the dose values.')

        # TODO: add more options for inputting the dose information
        form.addParam('inputDoseType',
                      params.EnumParam,
                      choices=['Scipion import', 'Fixed dose'],
                      default=SCIPION_IMPORT,
                      label='Input dose source',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Where to find the dose information:\n'
                           '- Scipion import: use the dose provided '
                           'during import of the tilt-series\n'
                           '- Fixed dose: manually input fixed dose '
                           'for each image of the input file, '
                           'in electrons/square Ångstrom.')

        form.addParam('fixedImageDose',
                      params.FloatParam,
                      default=FIXED_DOSE,
                      label='Fixed dose (e/sq Å)',
                      condition='inputDoseType == %i' % FIXED_DOSE,
                      help='Fixed dose for each image of the input file, '
                           'in electrons/square Ångstrom.')

        form.addParam('processOddEven',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Filter odd/even',
                      help='If True, the full tilt series and the associated odd/even tilt series will be processed. '
                           'The applied dose for the odd/even tilt series will be exactly the same.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.doseFilterStep, ts.getObjId())
            self._insertFunctionStep(self.createOutputStep, ts.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def doseFilterStep(self, tsObjId):
        """Apply the dose filter to every tilt series"""

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        firstItem = ts.getFirstItem()

        paramsMtffilter = {
            'input': firstItem.getFileName(),
            'output': os.path.join(extraPrefix, firstItem.parseFileName()),
            'pixsize': ts.getSamplingRate(),
            'voltage': ts.getAcquisition().getVoltage(),
        }

        argsMtffilter = "-input %(input)s " \
                        "-output %(output)s " \
                        "-PixelSize %(pixsize)f " \
                        "-Voltage %(voltage)d "

        if self.initialDose.get() != 0.0:
            argsMtffilter += f"-InitialDose {self.initialDose.get():f} "

        if self.inputDoseType.get() == SCIPION_IMPORT:
            outputDoseFilePath = os.path.join(tmpPrefix,
                                              firstItem.parseFileName(extension=".dose"))

            utils.generateDoseFileFromDoseTS(ts, outputDoseFilePath)

            paramsMtffilter.update({
                'typeOfDoseFile': 2,
                'doseWeightingFile': outputDoseFilePath,
            })

            argsMtffilter += "-TypeOfDoseFile %(typeOfDoseFile)d " \
                             "-DoseWeightingFile %(doseWeightingFile)s "

        elif self.inputDoseType.get() == FIXED_DOSE:
            paramsMtffilter.update({
                'fixedImageDose': self.fixedImageDose.get()
            })

            argsMtffilter += "-FixedImageDose %(fixedImageDose)f"

        Plugin.runImod(self, 'mtffilter', argsMtffilter % paramsMtffilter)

        if self.applyToOddEven(ts):
            oddFn = firstItem.getOdd().split('@')[1]
            paramsMtffilter['input'] = oddFn
            paramsMtffilter['output'] = os.path.join(extraPrefix, tsId+EXT_MRCS_TS_ODD_NAME)

            Plugin.runImod(self, 'mtffilter', argsMtffilter % paramsMtffilter)
            evenFn = firstItem.getEven().split('@')[1]
            paramsMtffilter['input'] = evenFn
            paramsMtffilter['output'] = os.path.join(extraPrefix, tsId+EXT_MRCS_TS_EVEN_NAME)
            Plugin.runImod(self, 'mtffilter', argsMtffilter % paramsMtffilter)

    def createOutputStep(self, tsObjId):
        """Generate output filtered tilt series"""

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        output = self.getOutputSetOfTiltSeries(self.inputSetOfTiltSeries.get())
        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)

        output.append(newTs)

        ih = ImageHandler()

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True, copyTM=True)
            newTi.setAcquisition(tiltImage.getAcquisition())
            if self.applyToOddEven(ts):
                locationOdd = index + 1, (os.path.join(extraPrefix, tsId+EXT_MRCS_TS_ODD_NAME))
                locationEven = index + 1, (os.path.join(extraPrefix, tsId+EXT_MRCS_TS_EVEN_NAME))
                newTi.setOddEven([ih.locationToXmipp(locationOdd), ih.locationToXmipp(locationEven)])
            else:
                newTi.setOddEven([])

            locationTi = index + 1, (os.path.join(extraPrefix,
                                                  tiltImage.parseFileName()))
            newTi.setLocation(locationTi)
            newTs.append(newTi)
            newTs.update(newTi)

        newTs.write(properties=False)

        output.update(newTs)
        output.write()

        self._store()

    def closeOutputSetsStep(self):
        self.TiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.TiltSeries.write()
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        if self.inputDoseType.get() == SCIPION_IMPORT:
            for ts in self.inputSetOfTiltSeries.get():
                if ts.getFirstItem().getAcquisition().getDosePerFrame() is None:
                    validateMsgs.append("%s has no dose information stored "
                                        "in Scipion Metadata. To solve this, import "
                                        "tilt-series using the mdoc option." %
                                        ts.getTsId())

        return validateMsgs

    def _summary(self):
        summary = []

        summary.append("%d input tilt-series" % self.inputSetOfTiltSeries.get().getSize())

        if self.TiltSeries:
            summary.append("%d tilt-series dose-weighted" % self.TiltSeries.getSize())

        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("The dose-weighting has been applied to %d "
                           "tilt-series using the IMOD *mtffilter* command."
                           % self.TiltSeries.getSize())
        return methods
