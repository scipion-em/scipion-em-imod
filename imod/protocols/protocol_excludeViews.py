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
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.object import Set
import tomo.objects as tomoObj
from imod import Plugin
from imod.protocols.protocol_base import ProtImodBase, OUTPUT_TILTSERIES_NAME
from pyworkflow.protocol import LEVEL_ADVANCED


class ProtImodExcludeViews(ProtImodBase):
    """
    excludeviews - Reversibly remove views from a tilt series stack
    If you use this protocol, make sure tis output tilt series is use for everything else
    CTF estimation, per particle per tilt, tomogram reconstruction....
    More info:
        https://bio3d.colorado.edu/imod/doc/man/excludeviews.html
    """

    _label = 'Exclude views'
    _devStatus = BETA
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: tomoObj.SetOfTiltSeries}

    excludeViewsInfoMatrix = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('excludeViewsFile',
                      params.FileParam,
                      expertLevel=LEVEL_ADVANCED,
                      label='Exclude views file',
                      help='File containing the views to be excluded for each tilt-series belonging to the set.\n\n'
                           'The format of the text file must be two columns, the first one being the tilt series ID of '
                           'the series from which the views will be excluded and the second the views to exclude, '
                           'numbered from 1. The syntax for this exclude list is a comma separated list of ranges with '
                           'no spaces between them (e.g., 1,4-5,60-70). \n\n'
                           'An example of this file comes as follows:\n'
                           'stack1 1,4-6,8,44-47\n'
                           'stack2 3,10-12,24\n'
                           '...')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):

        for tsId in self.inputSetOfTiltSeries.get().getIdSet():
            self._insertFunctionStep(self.excludeViewsStep, tsId)
            self._insertFunctionStep(self.generateOutputStackStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ----------------------------
    def getExcludedViewsFromFile(self):
        if self.excludeViewsInfoMatrix is None:

            self.excludeViewsInfoMatrix = utils.readExcludeViewsFile(self.excludeViewsFile.get())

        return self.excludeViewsInfoMatrix

    def getExcludedViewsFromMatrix(self, ts):
        """ Returns the indexes of the tilt to exclude for a specified tilt series read from de input file """
        matrix = self.getExcludedViewsFromFile()

        pattern = matrix.get(ts.getTsId(), [])

        views = self.makeExclusionPatternAsList(pattern)

        # Return the matrix or an empty list
        return views

    def getExcludedViews(self, ts):
        """ Returns the indexes of the tilt to exclude for a specified tilt series"""
        if self.excludeViewsFile.get():
            self.info("Using excluded views from %s" % self.excludeViewsFile.get() )
            return self.getExcludedViewsFromMatrix(ts)
        else:
            return ts._getExcludedViewsIndex()


    def excludeViewsStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        # Make folder to work with the tilt series
        extraPrefix = self._getExtraPath(tsId)
        path.makePath(extraPrefix)

        # Get the ts stack filename and coy it
        firstItem = ts.getFirstItem()
        outputFileName = os.path.join(extraPrefix, firstItem.parseFileName())
        excludedViews = self.getExcludedViews(ts)

        if excludedViews:

            self.info("Excluded views (%s) found for %s." % (excludedViews, tsId))

            path.copyFile(firstItem.getFileName(), outputFileName)

            paramsAlignment = {
                'stackName': outputFileName,
                'viewsToExclude': ",".join([str(index) for index in excludedViews]),
            }

            argsAlignment = "-StackName %(stackName)s " \
                            "-ViewsToExclude %(viewsToExclude)s "

            Plugin.runImod(self, 'excludeviews', argsAlignment % paramsAlignment)
        else:
            # Just create the link
            self.info("No views to exclude for %s." % tsId)
            path.createLink(firstItem.getFileName(), outputFileName)

    def generateOutputStackStep(self, tsObjId):
        output = self.getOutputSetOfTiltSeries(self.inputSetOfTiltSeries.get())

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        output.append(newTs)

        excludedViews = self.getExcludedViews(ts)

        for index, tiltImage in enumerate(ts):
            stackIndex = index+1
            if not stackIndex in excludedViews:
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True, copyTM=True)
                newTi.setAcquisition(tiltImage.getAcquisition())
                newTi.setLocation(stackIndex, (os.path.join(extraPrefix, tiltImage.parseFileName())))
                newTs.append(newTi)
            else:
                self.info("%s excluded from %s." % (stackIndex, tsId))
        newTs.write(properties=False)

        output.update(newTs)
        output.updateDim()
        output.write()
        self._store()

    def closeOutputSetsStep(self):
        self.TiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.TiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------

    def makeExclusionPatternAsList(self, excludedViews):
        excludedViewsAsList = []

        vector = excludedViews.split(',')

        for element in vector:
            elementVector = element.split('-')

            if len(elementVector) > 1:
                for i in range(int(elementVector[0]), int(elementVector[1])  + 1):
                    excludedViewsAsList.append(int(i))
            else:
                excludedViewsAsList.append(int(elementVector[0]))

        return excludedViewsAsList

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if self.TiltSeries:
            summary.append("Excluded views:\n")

            for tsIn, tsOut in zip(self.inputSetOfTiltSeries.get(), self.TiltSeries):
                summary.append("Tilt-series ID: %s.   Size: %d ----> %d."
                               % (tsIn.getTsId(),
                                  tsIn.getSize(),
                                  tsOut.getSize()))
        else:
            summary.append("Output classes not ready yet.")

        return summary

    def _methods(self):
        methods = []
        if self.TiltSeries:
            methods.append("Excluded views:\n")

            for tsIn, tsOut in zip(self.inputSetOfTiltSeries.get(), self.TiltSeries):
                methods.append("Tilt-series ID: %s.   Size: %d ----> %d."
                               % (tsIn.getTsId(),
                                  tsIn.getSize(),
                                  tsOut.getSize()))
        else:
            methods.append("Output classes not ready yet.")

        return methods
