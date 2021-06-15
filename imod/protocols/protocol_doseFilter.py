# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# **************************************************************************

import os
import imod.utils as utils
import pwem.objects as data
from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from imod import Plugin
from pwem.emlib.image import ImageHandler


class ProtImodDoseFilter(EMProtocol, ProtTomoBase):
    """
    Tilt-series' dose filtering based on the IMOD procedure.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/mtffilter.html
    """

    _label = 'dose filter'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series to be filtered.')

        form.addParam('initialDose',
                      params.FloatParam,
                      default=0.0,
                      label='Initial dose (e/sq A)',
                      help='Dose applied before any of the images in the input file were taken; this value will be '
                           'added to all the prior dose values, however they were obtained.')

        form.addParam('useFixedDose',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Use fixed dose',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Use a fixed dose fo every image in the tilt series instead the obtained form the '
                           'acquisition')

        form.addParam('fixedImageDose',
                      params.FloatParam,
                      default=1.0,
                      label='Fixes dose (e/sq A)',
                      condition='useFixedDose==0',
                      help='Fixed dose for each image of the input file, in electrons/square Angstrom.')


    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.computeXcorrStep, ts.getObjId())
            self._insertFunctionStep(self.createOutputStep, ts.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)


    # --------------------------- STEPS functions ----------------------------
    def computeXcorrStep(self, tsObjId):
        """Apply the dose fitler to every tilt series"""

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        paramsMtffilter = {
            'input': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            'output': os.path.join(extraPrefix, ts.getFirstItem().parseFileName()),
        }

        argsMtffilter = "-input %(input)s " \
                        "-output %(output)s " \

        if self.useFixedDose.get() == 1:
            paramsMtffilter.update({
                'fixedImageDose': self.minimumViewsPhaseShift.get()
            })

            argsMtffilter += "-FixedImageDose %(fixedImageDose)f"

        Plugin.runImod(self, 'tiltxcorr', argsMtffilter % paramsMtffilter)


    def createOutputStep(self, tsObjId):
        """Generate output filtered tilt series"""

        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)

        outputSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1, (os.path.join(extraPrefix, tiltImage.parseFileName())))

            newTs.append(newTi)

        newTs.write(properties=False)

        outputSetOfTiltSeries.update(newTs)
        outputSetOfTiltSeries.write()

        self._store()

    def closeOutputSetsStep(self):
        self.getOutputSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfTiltSeries(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries(suffix="Filtered")
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfTiltSeries.getSize()))
        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           "Interpolated Tilt-Series: %d.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           % (self.outputSetOfTiltSeries.getSize()))
        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           "Also, interpolation has been completed for %d Tilt-series.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
