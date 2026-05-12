# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
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
import logging
import traceback
from typing import Union

import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from imod.protocols.protocol_base import IN_TS_SET
from pwem.convert.headers import setMRCSamplingRate
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr, redStr
from pyworkflow.utils.retry_streaming import retry_on_sqlite_lock
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from imod.protocols import ProtImodBase
from imod.constants import OUTPUT_TILTSERIES_NAME, ODD, EVEN, EXCLUDE_VIEWS_PROGRAM

logger = logging.getLogger(__name__)


class ProtImodExcludeViews(ProtImodBase):
    """
    Removes selected tilt images from a tilt-series in a reversible and
    workflow-safe manner using IMOD excludeviews procedures. The protocol
    is designed to generate cleaned tilt-series datasets that exclude
    unusable, corrupted, low-quality, or experimentally problematic views
    before downstream tomographic processing steps such as CTF estimation,
    particle extraction, subtomogram analysis, or tomogram reconstruction.
    More info:
        https://bio3d.colorado.edu/imod/doc/man/excludeviews.html

    AI Generated:

    Exclude Views (ProtImodExcludeViews) - User Manual
        Overview

        The Exclude Views protocol removes selected images from a tilt
        series while preserving the integrity and acquisition metadata of
        the remaining dataset. In cryo-electron tomography workflows, some
        tilt images may become unsuitable for reconstruction because of
        excessive drift, charging, contamination, ice thickness changes,
        beam-induced damage, failed alignment, missing signal, or detector
        artifacts. Excluding these problematic views improves the overall
        quality and stability of downstream processing.

        From a biological perspective, the protocol is particularly useful
        when individual tilts compromise reconstruction quality or produce
        artifacts that interfere with interpretation of cellular or
        macromolecular structures. Removing these images can significantly
        improve tomogram consistency, especially in high-tilt regions where
        signal-to-noise ratios are naturally lower.

        Inputs and General Workflow

        The protocol requires an input set of tilt series in which some
        views may already be marked as disabled. These disabled images are
        interpreted as candidates for exclusion. The protocol then creates
        a new tilt-series stack containing only the retained views while
        preserving the geometric and acquisition consistency of the dataset.

        In many practical cryo-ET workflows, excluded views are identified
        during earlier preprocessing or quality-control stages. For
        example, users may manually disable views after visual inspection,
        remove images with severe motion or charging, or exclude tilts with
        failed CTF estimation. The protocol consolidates these decisions
        into a clean and reconstruction-ready dataset.

        The resulting tilt series should normally replace the original one
        in all subsequent processing stages. Continuing reconstruction or
        alignment using the original stack after excluding views may lead
        to inconsistencies in geometry, dose weighting, or metadata
        interpretation.

        Biological and Experimental Considerations

        Excluding views inevitably reduces angular sampling. While removing
        poor-quality tilts often improves reconstruction quality overall,
        excessive exclusion can increase the missing wedge effect and reduce
        isotropy in the final tomogram. Users should therefore balance data
        quality against angular completeness.

        In practice, removing a small number of severely corrupted views is
        usually beneficial. However, excluding many neighboring tilts,
        particularly around critical angular ranges, may introduce
        directional artifacts or reduce structural interpretability. This
        consideration becomes especially important for subtomogram averaging
        and quantitative structural analysis.

        The protocol automatically updates the effective angular range of
        the dataset after exclusion. This ensures that downstream software
        correctly interprets the minimum and maximum tilt angles, total dose
        accumulation, and related acquisition parameters.

        Odd and Even Tilt-Series Handling

        The protocol can also operate on odd and even tilt-series datasets
        when these are available. This capability is useful in workflows
        involving independent half-set processing, denoising validation, or
        resolution assessment strategies. Maintaining synchronization
        between the excluded views of the full, odd, and even datasets helps
        preserve consistency during subsequent reconstruction and analysis.

        Metadata Consistency

        A key objective of the protocol is preserving acquisition coherence
        after view exclusion. The resulting tilt series maintains updated
        information about angular coverage, dose accumulation, and image
        ordering so that downstream algorithms receive physically meaningful
        acquisition metadata.

        This is particularly important for reconstruction software that
        relies on accurate tilt geometry and dose information. Incorrect or
        outdated metadata after removing views may lead to reconstruction
        artifacts, weighting errors, or inaccurate interpretation of sample
        geometry.

        Outputs and Their Interpretation

        After execution, the protocol produces a new set of tilt series in
        which excluded images have been removed from the stacks. The output
        datasets preserve the acquisition structure of the remaining views
        while reflecting the updated angular and dose characteristics of
        the experiment.

        If no views are excluded, the protocol preserves the original tilt
        series while maintaining compatibility with the rest of the
        processing workflow. This allows users to integrate the protocol
        into automated pipelines without introducing unnecessary changes
        when all views are acceptable.

        Practical Recommendations

        In routine cryo-ET processing, it is advisable to inspect tilt
        images carefully before reconstruction and disable only those views
        that clearly compromise data quality. Common indicators include
        severe drift, extreme defocus instability, contamination, charging,
        excessive beam damage, or strong alignment residuals.

        Users should avoid excluding views solely because they appear noisy,
        particularly at high tilts where lower contrast is expected. In many
        cases, retaining a noisy but geometrically correct tilt contributes
        more useful information than removing it.

        After exclusion, it is recommended to verify the angular coverage
        and visually inspect reconstructed tomograms for directional
        artifacts or missing information. When many views are removed,
        additional caution is needed during biological interpretation.

        Final Perspective

        Excluding problematic views is often a necessary quality-control
        step in cryo-electron tomography. Although the process reduces the
        total amount of data, careful removal of corrupted images can
        substantially improve tomogram quality, alignment stability, and
        downstream structural interpretation. The most reliable results are
        obtained when exclusion decisions are guided by both image quality
        and biological relevance rather than by automated filtering alone.
    """

    _label = 'Exclude views'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsDict = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(Message.LABEL_INPUT)
        super().addInTsSetFormParam(form)
        self.addOddEvenParams(form)
        form.addParallelSection(threads=2, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId in self.tsDict.keys():
            exclStepId = self._insertFunctionStep(self.excludeViewsStep,
                                                  tsId,
                                                  prerequisites=[],
                                                  needsGPU=False)
            outStepId = self._insertFunctionStep(self.createOutputStep,
                                                 tsId,
                                                 prerequisites=exclStepId,
                                                 needsGPU=False)
            closeSetStepDeps.append(outStepId)
        self._insertFunctionStep(self.closeOutputSetsStep,
                                 OUTPUT_TILTSERIES_NAME,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        super()._initialize()
        self.tsDict = {ts.getTsId(): ts.clone() for ts in self.getInputTsSet().iterItems()}

    def excludeViewsStep(self, tsId: str):
        self.genTsPaths(tsId)
        try:
            ts = self.tsDict[tsId]
            firstTi = ts.getFirstEnabledItem()
            tsFileName = firstTi.getFileName()
            outTsFn, outTsOddFn, outTsEvenFn = self.getTmpFileNames(tsId)
            self._runExcludeViews(ts, tsFileName, outTsFn)
            if self.doOddEven:
                # ODD
                logger.info(cyanStr(f'tsId = {tsId} {ODD}:'))
                tsFnOdd = firstTi.getOdd()
                self._runExcludeViews(ts, tsFnOdd, outTsOddFn, suffix=ODD)
                # EVEN
                logger.info(cyanStr(f'tsId = {tsId} {EVEN}:'))
                tsFnEven = firstTi.getEven()
                self._runExcludeViews(ts, tsFnEven, outTsEvenFn, suffix=EVEN)

        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> {EXCLUDE_VIEWS_PROGRAM} execution '
                                f'failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        ts = self.tsDict[tsId]
        if tsId in self.failedItems:
            self.addToOutFailedSet(ts)
            return

        try:
            outFn = self.getExtraOutFile(tsId)
            if ts.hasExcludedViews():  # Then a new binary will be generated. If not, the original file remains linked
                setMRCSamplingRate(outFn, ts.getSamplingRate())  # Update the apix value in file header
            self._registerOutput(ts, outFn)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    @retry_on_sqlite_lock(log=logger)
    def _registerOutput(self, ts: TiltSeries, outFn: str):
        tsId = ts.getTsId()
        with self._lock:
            outTsSet = self.getOutputSetOfTS(self.getInputTsSet(pointer=True))
            outTs = TiltSeries()
            outTs.copyInfo(ts)
            outTsSet.append(outTs)

            angleMin = 999
            angleMax = -999
            accumDose = 0
            initialDose = 999
            tiList = []
            for ti in ts.iterItems():
                if ti.isEnabled():
                    # Update the acquisition of the TS. The accumDose, angle min and angle max for the re-stacked TS, as
                    # these values may change if the removed tilt-images are the first or the last, for example.
                    tiAngle = ti.getTiltAngle()
                    angleMin = min(tiAngle, angleMin)
                    angleMax = max(tiAngle, angleMax)
                    accumDose = max(ti.getAcquisition().getAccumDose(), accumDose)
                    initialDose = min(ti.getAcquisition().getDoseInitial(), initialDose)

                    outTi = TiltImage()
                    outTi.copyInfo(ti)
                    outTi.setFileName(outFn)
                    self.setTsOddEven(tsId, outTi, binGenerated=True)
                    tiList.append(outTi)

            if ts.hasExcludedViews():
                # Update the acquisition minAngle and maxAngle values of the tilt-series
                acq = outTs.getAcquisition()
                acq.setAngleMin(angleMin)
                acq.setAngleMax(angleMax)
                acq.setAccumDose(accumDose)
                acq.setDoseInitial(initialDose)
                outTs.setAcquisition(acq)
                # Update the acquisition minAngle and maxAngle values of each tilt-image acq while preserving their
                # specific accum and initial dose values
                for tiOut in tiList:
                    tiAcq = tiOut.getAcquisition()
                    tiAcq.setAngleMin(angleMin)
                    tiAcq.setAngleMax(angleMax)
                    tiOut.setAcquisition(tiAcq)
                    outTs.append(tiOut)
                outTs.setAnglesCount(len(outTs))
            else:
                for tiOut in tiList:
                    outTs.append(tiOut)

            outTs.write()
            outTsSet.update(outTs)
            outTsSet.write()
            self._store(outTsSet)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        output = getattr(self, OUTPUT_TILTSERIES_NAME, None)
        if output is not None:
            summary.append("Excluded views:\n")

            tsInSet = self.getInputTsSet().iterItems(orderBy='_tsId')
            tsOutSet = output.iterItems(orderBy='_tsId')
            for tsIn, tsOut in zip(tsInSet, tsOutSet):
                summary.append(f"Tilt-series: {tsIn.getTsId()}; "
                               f"Size: {tsIn.getSize()} ---> {tsOut.getSize()}")
        else:
            summary.append("Outputs are not ready yet.")

        return summary

    # --------------------------- UTILS functions -----------------------------
    def _runExcludeViews(self,
                         ts: TiltSeries,
                         inTsFn: str,
                         outTsFn: str,
                         suffix: str = "") -> None:
        tsId = ts.getTsId()
        finalLocation = self.getExtraOutFile(tsId, suffix=suffix)
        if ts.hasExcludedViews():
            logger.info(cyanStr(f'tsId = {tsId} {suffix}: excluding the disabled views...'))
            # The original file is moved to tmp
            path.copyFile(inTsFn, outTsFn)
            excludeViewsInd = ts.getTsExcludedViewsIndices(ts.getTsPresentAcqOrders())
            eVParams = {
                '-StackName': outTsFn,
                '-ViewsToExclude': ",".join(map(str, excludeViewsInd)),
            }
            self.runProgram(EXCLUDE_VIEWS_PROGRAM, eVParams)
            # The generated file overrides the original one in tmp, and it is moved to extra
            path.moveFile(outTsFn, finalLocation)
        else:
            # Just create the link
            logger.info(cyanStr(f"tsId = {tsId} -> No views to exclude."))
            path.createLink(inTsFn, finalLocation)

