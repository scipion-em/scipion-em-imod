# *****************************************************************************
# *
# * Authors:     Scipion Team [1]
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
from typing import List
from imod.protocols import ProtImodBase
from imod.protocols.protocol_base import BINNING_FACTOR
from pyworkflow.protocol import params, Form

# Float densities modes
FLOAT_DENSITIES_CHOICES = ['No adjust',
                           'adjust each section to fill the data range',
                           'scaled to common mean and standard deviation',
                           'shifted to a common mean without scaling',
                           'Rescaled to a min and max']
# Float densities values
NO_ADJUST = 0
EACH_SECTION_FILL_RANGE = 1
SCALED_COMMON_STD_MEAN = 2
SHIFTED_COMMON_MEAN = 3
SHIFTED_TO_MINMAX = 4


class ProtImodBasePreprocess(ProtImodBase):
    """
    Preprocesses tilt-series images or tomograms before downstream
    reconstruction, alignment, or analysis workflows. The protocol
    provides a standardized framework for reducing dataset size,
    normalizing image intensities, adjusting storage precision, and
    improving visual consistency across datasets.

    AI Generated:

    IMOD Preprocessing Base (ProtImodBasePreprocess) - User Manual
        Overview

        The preprocessing protocol prepares tilt-series data or tomograms
        for subsequent cryo-electron tomography operations by applying
        controlled image normalization and resampling procedures. Its
        objective is to improve data consistency, reduce computational
        demands, and ensure compatibility with later alignment or
        reconstruction stages within tomography workflows.

        In practical cryo-ET processing, preprocessing is often one of
        the first biologically important stages because it influences
        how effectively fiducials, structural features, and membrane
        boundaries can be identified in later analyses. Appropriate
        preprocessing can substantially improve reconstruction quality,
        whereas unsuitable normalization or excessive reduction in image
        detail may compromise downstream interpretation.

        Binning and Dataset Reduction

        One of the main functions of the protocol is image binning.
        Binning reduces the dimensions of the input images or tomograms
        by averaging neighboring pixels together. This operation lowers
        memory usage and accelerates computational processing, making
        large tomography datasets easier to handle.

        From a biological perspective, binning is commonly used during
        exploratory stages, coarse alignment, or rapid screening of
        datasets. Lower-resolution data can often provide sufficient
        information for identifying large-scale structures, fiducial
        trajectories, or general sample organization. However, excessive
        binning removes high-resolution detail and may obscure small
        macromolecular complexes, thin membranes, or subtle conformational
        features.

        In many workflows, users initially process heavily binned data
        to optimize alignment parameters and later repeat critical steps
        using minimally binned or unbinned datasets for final analysis.

        Density Normalization Strategies

        The protocol supports several approaches for adjusting image
        densities and intensity distributions. These normalization
        procedures help produce more homogeneous datasets and improve
        stability in later correlation-based processing steps.

        A commonly recommended strategy is normalization to a shared
        mean and standard deviation across all sections. This approach
        creates consistent intensity statistics throughout the tilt
        series and is particularly useful when acquisition conditions
        vary between projections or when illumination changes across
        the dataset.

        Other normalization approaches focus on shifting mean intensity
        values or remapping density ranges. These methods may be useful
        for visualization or compatibility with external software tools,
        but they should be selected carefully because excessive intensity
        manipulation can distort biologically meaningful contrast.

        For cryo-ET applications involving quantitative interpretation,
        users should generally prefer conservative normalization methods
        that preserve relative signal variations while reducing strong
        acquisition artifacts.

        Storage Precision and Output Modes

        The protocol allows the selection of different output storage
        formats and numerical precision levels. This flexibility is
        important because tomography datasets are often extremely large,
        and storage optimization can significantly reduce disk usage and
        transfer times.

        Lower precision formats may be sufficient for rapid visualization
        or intermediate alignment procedures, whereas floating-point
        representations are generally preferable for quantitative image
        processing and high-quality reconstruction workflows.

        Biological users should recognize that aggressive reduction of
        precision can introduce rounding artifacts or reduce weak density
        information, particularly in low-contrast cryo-ET datasets.

        Antialias Filtering and Image Quality

        During image reduction, antialias filtering helps preserve image
        quality by minimizing artifacts introduced by resampling. The
        protocol provides several filtering strategies with different
        balances between smoothing, sharpness preservation, and artifact
        suppression.

        In biological tomography, preserving structural continuity is
        particularly important for membranes, filament systems, and
        macromolecular assemblies with fine edges. Strong smoothing may
        reduce noise but can also blur biologically meaningful features.
        Conversely, sharper filters preserve detail but may introduce
        ringing artifacts near high-contrast boundaries.

        The optimal antialiasing strategy depends on specimen type,
        imaging conditions, and the intended downstream application.
        Exploratory workflows often tolerate stronger smoothing, whereas
        subtomogram averaging or segmentation pipelines may benefit from
        more detail-preserving filters.

        Validation and Recommended Usage

        The protocol includes validation safeguards to discourage the use
        of outdated or biologically suboptimal normalization approaches.
        This helps maintain compatibility with modern cryo-ET processing
        standards and encourages more stable preprocessing practices.

        For most routine workflows, moderate binning combined with
        normalization to a common mean and standard deviation provides
        a robust starting point. Users should visually inspect processed
        datasets to confirm that membrane continuity, fiducial contrast,
        and structural integrity remain biologically meaningful after
        preprocessing.

        Final Perspective

        Preprocessing is not merely a technical preparation stage but an
        important determinant of downstream tomographic quality. Careful
        control of binning, normalization, filtering, and storage
        precision can substantially improve the reliability of alignment,
        reconstruction, segmentation, and interpretation in cryo-electron
        tomography workflows.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _defineParams(self,
                      form: Form,
                      isTomogram: bool = False):
        objStr = 'tomograms' if isTomogram else 'images'

        form.addParam(BINNING_FACTOR,
                      params.IntParam,
                      default=1,
                      label='Binning',
                      important=True,
                      validators=[params.GE(1)],
                      help=f'Binning is an scaling factor for the output {objStr}. It must be an integer greater than '
                           f'1. IMOD uses ordinary binning to reduce {objStr} in size by the given factor. The value '
                           f'of a binned pixel is the average of pixel values in each block of pixels being binned. '
                           f'Binning is applied before all.')

        form.addParam('floatDensities',
                      params.EnumParam,
                      choices=FLOAT_DENSITIES_CHOICES,
                      default=2,
                      label='Adjust densities mode',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help=f'Adjust densities of sections individually. Modes:\n\n'
                           f'*0 - No adjustment performed*\n\n'

                           f'*1 - Adjust each section to fill the data range*\n\n'

                      # '-Range between min and max: This option will scale the gray values'
                      # 'to be in a range given by a minimum and a maximum values.'
                      # 'This is the mode 1 in newstack flag -floatDensities.\n'

                           f'*2 - Scaled to common mean and standard deviation*:\n\n'
                           f'This is the most '
                           f'common normalization procedure. The new tilt series will have'
                           f'a meand and a standard deviation introduced by the user. Generaly,'
                           f'a zero mean and a standard deviation one is a good choice.'
                           f'This is the mode 2 in newstack flag -floatDensities.\n\n'

                           f'*3 - Shifted to a common mean without scaling*:\n\n'
                           f'This option only'
                           f'add an offset to the gray values of the {objStr}. The offset will'
                           f'be calculated such as the new {objStr} will present a mean gray value'
                           f'introduced by the user. This is the mode 3 in newstack flag '
                           f'floatDensities.\n\n'

                           f'*4 - shifted to mean and rescaled to a min and max*:\n\nIn this case, an '
                           f'offset is added to the {objStr} in order to achieve a mean gray value '
                           f'then they are rescale the resulting minimum and maximum densities '
                           f'to the Min and Max values specified. This is the mode 4 in newstack '
                           f'flag -floatDensities.')

        # NEWSTACK - The -scale, -contrast, -multadd, and -float options are mutually exclusive except with -float 4
        # scaleRangeToggleCond = "floatDensities in [0, 4]"
        # -meansd (-mea) OR -MeanAndStandardDeviation   Two floats
        #               Scale all images to the given mean and standard deviation.  This
        #               option implies -float 2 and is incompatible with all other scaling options.
        floatDensMode2Cond = ('floatDensities in [%i, %i]' %
                              (SCALED_COMMON_STD_MEAN, SHIFTED_COMMON_MEAN))
        groupMeanSd = form.addGroup('Mean and SD',
                                    condition=floatDensMode2Cond,
                                    help=f'Scale all {objStr} to the given mean and standard deviation.')
        groupMeanSd.addParam('scaleMean',
                             params.FloatParam,
                             condition=floatDensMode2Cond,
                             default=0,
                             label='Mean',
                             help='Mean value for the rescaling')

        groupMeanSd.addParam('scaleSd',
                             params.FloatParam,
                             condition='floatDensities == %i ' % SCALED_COMMON_STD_MEAN,
                             default=1,
                             label='SD',
                             help='Standard deviation value for the rescaling')

        groupScale = form.addGroup('Scaling values',
                                   condition='floatDensities in [%i]' % (SHIFTED_TO_MINMAX))
        msg = 'This option will rescale the densities of all sections by the same factors so that the original ' \
              'minimum and maximum density will be mapped to the Min and Max values that are entered.'
        groupScale.addParam('scaleMax',
                            params.FloatParam,
                            default=255.,
                            label='Max.',
                            help=f'Maximum value for the rescaling. {msg}')

        groupScale.addParam('scaleMin',
                            params.FloatParam,
                            default=0.,
                            label='Min.',
                            help=f'Minimum value for the rescaling. {msg}')

        form.addParam('modeToOutput',
                      params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      choices=['default', '4-bit', '8-bit', 'signed 16-bit',
                               'unsigned 16-bit', '32-bit float'],
                      default=0,
                      label='Storage data type',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='The storage mode of the output file. The '
                           'default is the mode of the first input file, '
                           'except for a 4-bit input file, where the default '
                           'is to output as bytes')

        form.addParam('antialias',
                      params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      choices=['None', 'Blackman', 'Triangle', 'Mitchell',
                               'Lanczos 2 lobes', 'Lanczos 3 lobes'],
                      default=5,
                      label='Antialias method:',
                      display=params.EnumParam.DISPLAY_COMBO,
                      help=f'Type of antialiasing filter to use when reducing {objStr}.\n'
                           f'The available types of filters are:\n\n'
                           f'None - Antialias will not be applied\n'
                           f'Blackman - fast but not as good at antialiasing as slower filters\n'
                           f'Triangle - fast but smooths more than Blackman\n'
                           f'Mitchell - good at antialiasing, smooths a bit\n'
                           f'Lanczos 2 lobes - good at antialiasing, less smoothing than Mitchell\n'
                           f'Lanczos 3 lobes - slower, even less smoothing but more risk of ringing\n'
                           f'The default is Lanczos 3 as of IMOD 4.7. Although '
                           f'many people consider Lanczos 2 the best compromise '
                           f'among the various factors, that sentiment may be '
                           f'based on {objStr} of natural scenes where there are '
                           f'sharp edges.')

        self.addOddEvenParams(form, isTomogram=isTomogram)
        form.addParallelSection(threads=3, mpi=0)

    # --------------------------- INFO functions ------------------------------
    def _validate(self) -> List[str]:
        errorMsg = []
        floatDensitiesMode = self.floatDensities.get()
        if floatDensitiesMode == EACH_SECTION_FILL_RANGE:
            errorMsg.append('Float densities mode "Adjust each section to fill the data range" has been '
                            'deprecated. Please choose another method.')

        if floatDensitiesMode == SHIFTED_COMMON_MEAN:
            errorMsg.append('Float densities mode "Shifted to a common mean without scaling" has been '
                            'deprecated. Please choose another method.')
        return errorMsg
