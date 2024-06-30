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

import pyworkflow.protocol.params as params
from pwem.emlib.image import ImageHandler as ih
from tomo.objects import TiltSeries, TiltImage, SetOfTiltSeries

from imod import utils
from imod.protocols.protocol_base import ProtImodBase
from imod.constants import XF_EXT, ODD, EVEN, OUTPUT_TILTSERIES_NAME


class ProtImodApplyTransformationMatrix(ProtImodBase):
    """
    Compute the interpolated tilt-series from its transform matrix.
    The protocol makes use of the IMod command newstack
    More info:
        https://bio3d.colorado.edu/imod/doc/man/newstack.html

    Generally, the tilt series has an associated transformation matrix
    which contains the alignment information. The transformation matrix
    is usually associated but not applied to avoid to accumulate interpolation
    errors during the image processing. This protocol allows to apply
    the transformation matrix to the tilt series
    """

    _label = 'Apply transformation'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series to apply the transformation matrix')

        form.addParam('binning', params.IntParam,
                      default=1,
                      label='Binning for the interpolated',
                      help='Binning to be applied to the interpolated tilt-series '
                           'in IMOD convention.\nBinning is a scaling factor '
                           'given by an integer greater than 1. IMOD uses ordinary '
                           'binning (with antialiasing filter) to reduce images in '
                           'size by the given factor. The value of a binned pixel '
                           'is the average of pixel values in each block of pixels '
                           'being binned. Binning is applied before all other image '
                           'transformations.')

        form.addParam('taperInside',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Taper inwards from the edge?',
                      help='When the image is transformed areas with no information '
                           'are filled in (e.g. because of rotation).'
                           'Decide whether tapering is done inwards or outwards '
                           'from the edge.')
        
        form.addParam('linear',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=False,
                      label='Linear interpolation?',
                      help='From newstack man page: Use linear instead of cubic '
                           'interpolation to transform images. Linear interpolation '
                           'is more suitable when images are very noisy, but cubic '
                           'interpolation will preserve fine detail better when '
                           'noise is not an issue.')
        
        form.addParam('processOddEven',
                      params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      default=True,
                      label='Apply to odd/even',
                      help='If True, the full tilt series and the associated '
                           'odd/even tilt series will be processed. The '
                           'transformations applied to the odd/even tilt-series '
                           'will be exactly the same.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tsId in self.tsDict.keys():
            self._insertFunctionStep(self.computeAlignmentStep, tsId)
            self._insertFunctionStep(self.generateOutputStackStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ------------------------------
    def _initialize(self):
        self.tsDict = {ts.getTsId(): ts.clone() for ts in self.getInputSet()}

    def computeAlignmentStep(self, tsId):
        try:
            ts = self.tsDict[tsId]
            firstItem = ts.getFirstItem()
            self.genTsPaths(tsId)
            utils.genXfFile(ts, self.getExtraOutFile(tsId, ext=XF_EXT))

            params = self.getBasicNewstackParams(ts,
                                                 self.getExtraOutFile(tsId),
                                                 firstItem=firstItem,
                                                 xfFile=self.getExtraOutFile(tsId, ext=XF_EXT),
                                                 binning=self.binning.get(),
                                                 doSwap=True,
                                                 tsExcludedIndices=ts.getExcludedViewsIndex(),
                                                 doTaper=True)
            params["-taper"] = "1,1" if self.taperInside else "1,0"

            if self.linear:
                params["-linear"] = ""

            self.runProgram("newstack", params)

            if self.applyToOddEven(ts):
                oddFn = firstItem.getOdd().split('@')[1]
                evenFn = firstItem.getEven().split('@')[1]
                params['-input'] = oddFn
                params['-output'] = self.getExtraOutFile(tsId, suffix=ODD)
                self.runProgram("newstack", params)

                params['-input'] = evenFn
                params['-output'] = self.getExtraOutFile(tsId, suffix=EVEN)
                self.runProgram("newstack", params)

        except Exception as e:
            self._failedTs.append(tsId)
            self.error(f'Newstack execution failed for tsId {tsId} -> {e}')

    def generateOutputStackStep(self, tsId):
        ts = self.tsDict[tsId]
        if tsId in self._failedTs:
            self.createOutputFailedSet(ts)
        else:
            outputLocation = self.getExtraOutFile(tsId)
            if os.path.exists(outputLocation):
                output = self.getOutputInterpolatedTS(self.getInputSet())

                newTs = TiltSeries(tsId=tsId)
                newTs.copyInfo(ts)
                newTs.setInterpolated(True)
                acq = newTs.getAcquisition()
                acq.setTiltAxisAngle(0.)  # 0 because TS is aligned
                newTs.setAcquisition(acq)
                output.append(newTs)

                oddEvenFlag = self.applyToOddEven(ts)
                for index, tiltImage in enumerate(ts):
                    if tiltImage.isEnabled():
                        newTi = TiltImage()
                        newTi.copyInfo(tiltImage, copyId=False, copyTM=False)
                        acq = tiltImage.getAcquisition()
                        acq.setTiltAxisAngle(0.)
                        newTi.setAcquisition(acq)
                        newTi.setLocation(index+1, outputLocation)
                        if oddEvenFlag:
                            locationOdd = index+1, (self.getExtraOutFile(tsId, suffix=ODD))
                            locationEven = index+1, (self.getExtraOutFile(tsId, suffix=EVEN))
                            newTi.setOddEven([ih.locationToXmipp(locationOdd),
                                              ih.locationToXmipp(locationEven)])
                        else:
                            newTi.setOddEven([])

                        newTi.setSamplingRate(self._getOutputSampling())
                        newTs.append(newTi)

                dims = self._getOutputDim(outputLocation)
                newTs.setDim(dims)

                newTs.write(properties=False)
                output.update(newTs)
                output.write()
                self._store(output)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        for ts in self.getInputSet():
            if not ts.hasAlignment():
                validateMsgs.append("Some tilt-series from the input set "
                                    "are missing a transformation matrix.")
                break

        return validateMsgs

    def _summary(self):
        summary = []
        if self.InterpolatedTiltSeries:
            summary.append(f"Input tilt-series: {self.getInputSet().getSize()}\n"
                           f"Interpolations applied: {self.InterpolatedTiltSeries.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if self.InterpolatedTiltSeries:
            methods.append("The interpolation has been computed for "
                           f"{self.InterpolatedTiltSeries.getSize()} "
                           "tilt-series using the IMOD *newstack* command.")
        return methods

    # --------------------------- UTILS functions -----------------------------
    def _getOutputSampling(self) -> float:
        return self.getInputSet().getSamplingRate() * self.binning.get()
