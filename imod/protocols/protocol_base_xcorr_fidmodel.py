# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es) [1]
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
from pwem.protocols import EMProtocol
from pyworkflow.protocol import LEVEL_ADVANCED, NumericListParam, IntParam, LEVEL_NORMAL, FloatParam


class ProtImodBaseXcorrFidModel(EMProtocol):
    """
    Provides a common framework for configuring cross-correlation and fiducial-based
    alignment strategies in electron tomography workflows. The protocol focuses on
    defining biologically meaningful preprocessing and correlation settings that help
    improve the robustness of tilt-series alignment, especially in challenging datasets
    affected by noise, uneven contrast, large-scale gradients, or specimen deformation.

    AI Generated:

    Cross-Correlation Fiducial Model Base (ProtImodBaseXcorrFidModel) - User Manual
        Overview

        This protocol serves as a foundational configuration layer for tomography
        alignment procedures that rely on cross-correlation analysis and fiducial
        tracking approaches. Its purpose is to help users define the spatial regions
        and frequency filtering conditions that guide alignment algorithms toward the
        most reliable structural information within a tilt series.

        In cryo-electron tomography and related imaging workflows, alignment quality
        strongly determines the interpretability of the final reconstruction. Poorly
        aligned projections can produce blurred tomograms, distort structural details,
        and reduce the visibility of biologically meaningful features. This protocol
        therefore emphasizes preprocessing strategies that stabilize correlation-based
        alignment and improve fiducial localization.

        Trimming and Correlation Regions

        One of the most important concepts in tilt-series alignment is the selection
        of the image area used during correlation calculations. In many experimental
        datasets, image borders contain artifacts such as shadows, damaged regions,
        contamination, detector defects, or excessive noise. Excluding these areas
        often improves alignment stability and reduces the influence of non-biological
        signal.

        The protocol allows users to define trimming regions along the X and Y axes.
        These settings are especially useful for patch-tracking workflows, where local
        image regions are compared across tilts to estimate motion and geometric
        transformations. By restricting the analysis to informative areas of the field
        of view, users can obtain more stable correlations and reduce alignment errors.

        From a biological perspective, the selected regions should ideally contain
        recognizable structural information distributed consistently across the tilt
        series. Large empty solvent regions or highly damaged image areas should
        generally be avoided because they contribute little meaningful information to
        the alignment process.

        Frequency Filtering and Noise Management

        The protocol also provides configuration parameters for high-pass and low-pass
        frequency filtering. These filters are essential for controlling how image
        features contribute to cross-correlation calculations.

        High-pass filtering suppresses very large-scale intensity variations that may
        arise from uneven illumination, ice thickness gradients, contamination, or
        other low-frequency artifacts. In practical cryo-ET workflows, modest
        high-pass filtering often improves the ability of the alignment procedure to
        focus on structural details instead of broad intensity fluctuations.

        Low-pass filtering reduces the influence of high-frequency noise. This becomes
        particularly important in noisy datasets acquired under low-dose conditions,
        where detector noise and limited signal-to-noise ratios can destabilize
        correlation measurements. Proper low-pass filtering helps emphasize coherent
        structural signal while suppressing random fluctuations.

        The balance between high-pass and low-pass filtering is biologically important.
        Excessive filtering may remove genuine structural information, whereas
        insufficient filtering may allow artifacts or noise to dominate the alignment.
        Optimal values often depend on specimen size, acquisition conditions, detector
        quality, and binning level.

        Practical Use in Tomography Workflows

        In routine tomography processing, these preprocessing parameters are commonly
        adjusted during early optimization stages of a workflow. Users often begin
        with conservative filtering and moderate trimming values, then refine the
        parameters after inspecting alignment quality and reconstruction consistency.

        Patch-tracking approaches typically benefit from some degree of trimming
        because border regions frequently contain distortions that interfere with local
        correlation measurements. Fiducial-based workflows may also benefit from
        filtering adjustments when fiducial markers are difficult to distinguish from
        background noise.

        Datasets containing thick ice, crowded cellular environments, or highly noisy
        projections often require more careful tuning of the filtering parameters.
        Conversely, high-quality datasets with strong fiducial contrast may require
        only minimal preprocessing.

        Interpretation and Biological Considerations

        The parameters configured through this protocol do not directly alter the
        biological specimen itself, but they strongly influence how accurately the
        specimen geometry is reconstructed. Small differences in preprocessing can
        significantly affect the final tomogram quality, especially in subtomogram
        averaging or high-resolution structural studies.

        Users should therefore interpret these settings as part of the scientific
        optimization process rather than purely technical adjustments. The goal is to
        preserve biologically meaningful signal while minimizing artifacts that could
        bias alignment or reconstruction.

        Final Perspective

        Reliable tilt-series alignment depends not only on the alignment algorithm but
        also on the quality of the image regions and frequency information provided to
        that algorithm. Careful trimming of problematic image areas together with
        balanced frequency filtering can substantially improve alignment robustness,
        reconstruction clarity, and downstream biological interpretation in electron
        tomography workflows.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
    # -------------------------- DEFINE param functions -----------------------
    @staticmethod
    def addTrimingParams(form, pxTrimCondition=False, correlationCondition=True,
                         levelType=LEVEL_ADVANCED):
        """
        Generally, this form will be integrated in a groupForm,
        the group form argument is form. A set of flags
        control what elements are shown
        """
        form.addParam('pxTrim',
                      NumericListParam,
                      condition=pxTrimCondition,
                      label='Pixels to trim (x y without coma separator)',
                      default="40 40",
                      help='Pixels to trim off each side in X and Y.\n'
                           'Some trimming should be used for patch tracking',
                      expertLevel=levelType)

        xtrimming = form.addLine('Pixels to do correlation along X-axis',
                                 expertLevel=levelType,
                                 condition=correlationCondition,
                                 help="Starting and ending X coordinates of "
                                      "a region to correlate, based on the "
                                      "position of the region at zero tilt.")

        xtrimming.addParam('xmin',
                           IntParam,
                           label='X axis min (left)',
                           allowsNull=True,
                           expertLevel=levelType)

        xtrimming.addParam('xmax',
                           IntParam,
                           label='X axis max (right)',
                           allowsNull=True,
                           expertLevel=levelType)

        ytrimming = form.addLine('Pixels to do correlation along Y-axis',
                                 expertLevel=levelType,
                                 condition=correlationCondition,
                                 help="Starting and ending Y coordinates "
                                      "of a region to correlate, based on "
                                      "the position of the region at zero tilt.")

        ytrimming.addParam('ymin',
                           IntParam,
                           label='Y axis min (top)',
                           allowsNull=True,
                           expertLevel=levelType)

        ytrimming.addParam('ymax',
                           IntParam,
                           label='Y axis max (botton)',
                           allowsNull=True,
                           expertLevel=levelType)

    @staticmethod
    def filteringParametersForm(form, condition, levelType=LEVEL_NORMAL):
        filtering = form.addGroup('Filtering parameters',
                                  condition=condition,
                                  expertLevel=levelType)

        line1 = filtering.addLine('High pass filter',
                                  expertLevel=levelType,
                                  help="Some high pass filtering, using a small value of Sigma1 such "
                                       "as 0.03, may be needed to keep the program from being misled by very "
                                       "large scale features in the images.  If the images are noisy, some low "
                                       "pass filtering with Sigma2 and Radius2 is appropriate (e.g. 0.05 for "
                                       " Sigma2, 0.25 for Radius2).  If the images are binned, these values "
                                       "specify frequencies in the binned image, so a higher cutoff (less filtering) "
                                       "might be appropriate.\n\n"
                                       "*FilterRadius1*: Low spatial frequencies in the cross-correlation "
                                       "will be attenuated by a Gaussian curve that is 1 "
                                       "at this cutoff radius and falls off below this "
                                       "radius with a standard deviation specified by "
                                       "FilterSigma2. Spatial frequency units range from "
                                       "0 to 0.5.\n"
                                       "*Filter sigma 1*: Sigma value to filter low frequencies in the "
                                       "correlations with a curve that is an inverted "
                                       "Gaussian.  This filter is 0 at 0 frequency and "
                                       "decays up to 1 with the given sigma value. "
                                       "However, if a negative value of radius1 is entered, "
                                       "this filter will be zero from 0 to "
                                       "|radius1| then decay up to 1.")

        line1.addParam('filterRadius1',
                       FloatParam,
                       label='Filter radius 1',
                       default=0.0,
                       expertLevel=levelType)

        line1.addParam('filterSigma1',
                       FloatParam,
                       label='Filter sigma 1',
                       default=0.03,
                       expertLevel=levelType)

        line2 = filtering.addLine('Low pass filter',
                                  expertLevel=levelType,
                                  help="If the images are noisy, some low "
                                       "pass filtering with Sigma2 and Radius2 is appropriate (e.g. 0.05 for "
                                       " Sigma2, 0.25 for Radius2).  If the images are binned, these values "
                                       "specify frequencies in the binned image, so a higher cutoff (less filtering) "
                                       "might be appropriate.\n\n"
                                       "*Filter radius 2*: High spatial frequencies in the cross-correlation "
                                       "will be attenuated by a Gaussian curve that is 1 "
                                       "at this cutoff radius and falls off above this "
                                       "radius with a standard deviation specified by "
                                       "FilterSigma2.\n"
                                       "*Filter sigma 2*: Sigma value for the Gaussian rolloff below and "
                                       "above the cutoff frequencies specified by "
                                       "FilterRadius1 and FilterRadius2")

        line2.addParam('filterRadius2',
                       FloatParam,
                       label='Filter radius 2',
                       default=0.25,
                       expertLevel=levelType)

        line2.addParam('filterSigma2',
                       FloatParam,
                       label='Filter sigma 2',
                       default=0.05,
                       expertLevel=levelType)