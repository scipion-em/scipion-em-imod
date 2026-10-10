# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# **************************************************************************
from pyworkflow.protocol.constants import STEPS_SERIAL

from tomo.objects import SetOfCTFTomoSeries

from imod.protocols import ProtImodAutomaticCtfEstimation
from imod.constants import OUTPUT_CTF_SERIE


class ProtImodManualCtfEstimation(ProtImodAutomaticCtfEstimation):
    """
    Performs manual contrast transfer function estimation for tilt-series
    datasets using the IMOD ctfplotter workflow. The protocol provides an
    interactive environment in which users can visually inspect power spectra,
    evaluate defocus fitting quality, and manually save the resulting defocus
    estimations for downstream tomographic processing. More info:
    https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html

    AI Generated:

    Manual CTF Estimation (ProtImodManualCtfEstimation) - User Manual
        Overview

        The Manual CTF Estimation protocol provides an interactive framework
        for estimating contrast transfer function parameters from cryo-electron
        tomography tilt series using the IMOD ctfplotter utility. Its primary
        objective is to help users determine accurate defocus values through
        direct visual inspection of rotationally averaged power spectra.

        In cryo-electron microscopy, the contrast transfer function strongly
        influences image contrast and the preservation of structural
        information across spatial frequencies. Reliable CTF estimation is
        therefore essential for achieving accurate tomographic reconstruction,
        improving signal restoration, and supporting downstream analyses such
        as subtomogram averaging or particle detection.

        Unlike fully automated approaches, this protocol is specifically
        designed for workflows in which user supervision remains important.
        Interactive inspection is often beneficial for difficult datasets,
        especially when images contain low contrast, thick ice, strong
        astigmatism, uneven illumination, contamination, or complex biological
        backgrounds.

        Biological Motivation for Manual Inspection

        Tomographic tilt series frequently exhibit highly variable image
        quality across projection angles. At high tilts, increased sample
        thickness and reduced signal-to-noise ratios may complicate automated
        defocus fitting. In addition, biological specimens such as crowded
        cellular environments or membrane-rich regions may produce power
        spectra that are difficult to interpret automatically.

        Manual inspection allows users to determine whether the detected
        Thon rings are biologically meaningful and whether the fitted defocus
        values are consistent with the quality of the experimental data.
        Researchers can selectively include or exclude projection images during
        spectral averaging in order to improve fitting stability.

        This flexibility is particularly valuable when working with
        heterogeneous datasets, damaged projections, or tilt series affected by
        charging, drift, or radiation damage.

        Power Spectrum Averaging and Interpretation

        The protocol operates through rotationally averaged power spectra
        derived from overlapping image regions. Averaging multiple local
        spectra improves the visibility of oscillatory CTF features while
        reducing random noise contributions.

        From a biological perspective, the visibility of Thon rings provides
        indirect information about data quality and imaging conditions.
        Well-defined rings usually indicate stable imaging conditions and good
        preservation of high-resolution information, whereas weak or distorted
        rings may reflect specimen movement, charging effects, contamination,
        or limitations imposed by thick specimens.

        Users should interpret the fitted curves carefully, especially in
        tomographic datasets where anisotropic resolution and varying specimen
        thickness can complicate the appearance of the power spectrum.

        Interactive Workflow

        The protocol is designed around user interaction rather than automated
        execution. Researchers are expected to inspect spectra visually,
        evaluate fitting reliability, and manually confirm or save the final
        defocus estimations.

        This interactive philosophy is particularly useful for expert users who
        wish to maintain strict control over reconstruction quality. In many
        tomography workflows, small improvements in CTF estimation accuracy can
        significantly influence downstream alignment, reconstruction quality,
        and subtomogram interpretability.

        The protocol also allows users to process multiple tilt series within a
        unified interactive environment, supporting iterative inspection and
        quality control across complete experiments.

        Outputs and Downstream Applications

        The resulting output consists of CTF estimation information associated
        with the processed tilt series. These defocus parameters can later be
        used for CTF correction, phase flipping, dose weighting, or advanced
        tomographic refinement procedures.

        Accurate defocus estimation is especially important in workflows aiming
        for high-resolution subtomogram averaging, where errors in CTF modeling
        can limit the recovery of structural detail.

        Biological users should remember that the quality of downstream
        tomograms is closely linked to the reliability of the estimated CTF
        parameters. Careful manual validation is therefore often preferable to
        relying entirely on automated fitting in challenging datasets.

        Practical Recommendations

        For most cryo-electron tomography datasets, it is advisable to inspect
        several representative projections across the tilt range before
        accepting fitted values. High-tilt images commonly display weaker CTF
        oscillations and may require more cautious interpretation.

        Excluding projections with severe contamination, motion blur, or poor
        spectral quality often improves the stability of the averaged spectrum.
        Similarly, datasets collected with thick ice or crowded cellular
        material may benefit from careful selection of regions contributing to
        the spectral average.

        Researchers should also verify that estimated defocus trends remain
        biologically and experimentally consistent throughout the tilt series.

        Final Perspective

        Manual CTF estimation remains an important component of cryo-electron
        tomography workflows because many biological datasets still challenge
        fully automated procedures. Interactive validation allows researchers
        to combine computational fitting with expert interpretation of spectral
        quality and imaging conditions.

        The Manual CTF Estimation protocol provides a controlled environment
        for obtaining reliable defocus measurements while preserving the
        flexibility needed for complex and heterogeneous tomographic datasets.
    """

    _label = 'CTF estimation (manual)'
    _interactiveMode = True
    _possibleOutputs = {OUTPUT_CTF_SERIE: SetOfCTFTomoSeries}
    stepsExecutionMode = STEPS_SERIAL

    def __init__(self, **args):
        ProtImodAutomaticCtfEstimation.__init__(self, **args)
        self.OUTPUT_PREFIX = OUTPUT_CTF_SERIE

    @classmethod
    def worksInStreaming(cls):
        return False

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.runCTFEtimationStep, interactive=True)

    # --------------------------- STEPS functions ----------------------------

    def runCTFEtimationStep(self):
        from imod.viewers import ImodGenericView
        tsSet = self.getInputTsSet()
        view = ImodGenericView(None, self, tsSet,
                               createSetButton=True,
                               isInteractive=True)
        view.show()
        self.createOutput()

    def runAllSteps(self, obj):
        tsId = obj.getTsId()
        self.convertInputStep(tsId)
        expDefoci = self.getExpectedDefocus()
        self.ctfEstimation(tsId, expDefoci)

    def createOutput(self):
        suffix = self._getOutputSuffix(SetOfCTFTomoSeries)
        outputSetName = self.OUTPUT_PREFIX + str(suffix)
        for tsId in self.tsDict.keys():
            self.createOutputStep(tsId, outputSetName)

        self.closeOutputSetsStep()

    def _summary(self):
        summary = []
        return summary
