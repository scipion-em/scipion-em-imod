# *****************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Centro Nacional de Biotecnologia, CSIC, Spain
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


IMOD_HOME = 'IMOD_HOME'
IMOD_VIEWER_BINNING = 'IMOD_VIEWER_BINNING'
ETOMO_CMD = 'etomo'

VERSION_4_11_24 = '4.11.24'
VERSION_4_11_25 = '4.11.25'
VERSIONS = [VERSION_4_11_24, VERSION_4_11_25]
DEFAULT_VERSION = VERSION_4_11_25

# protocol constants below

SCIPION_IMPORT = 0
FIXED_DOSE = 1

# Fiducial model types
FIDUCIAL_MODEL = 0
PATCH_TRACKING = 1
# Patch tracking layout options
PT_FRACTIONAL_OVERLAP = 0
PT_NUM_PATCHES = 1

OUTPUT_TS_COORDINATES_NAME = "TiltSeriesCoordinates"
OUTPUT_FIDUCIAL_NO_GAPS_NAME = "FiducialModelNoGaps"
OUTPUT_FIDUCIAL_GAPS_NAME = "FiducialModelGaps"
OUTPUT_TILTSERIES_NAME = "TiltSeries"
OUTPUT_PREALI_TILTSERIES_NAME = "PreAlignedTiltSeries"
OUTPUT_ALI_TILTSERIES_NAME = "AlignedTiltSeries"
OUTPUT_TS_INTERPOLATED_NAME = "InterpolatedTiltSeries"
OUTPUT_TS_FAILED_NAME = "FailedTiltSeries"
OUTPUT_TOMOS_FAILED_NAME = "FailedTomograms"
OUTPUT_CTF_SERIE = "CTFTomoSeries"
OUTPUT_TOMOGRAMS_NAME = "Tomograms"
OUTPUT_COORDINATES_3D_NAME = "Coordinates3D"
EXT_MRCS_TS_EVEN_NAME = "even.mrcs"
EXT_MRCS_TS_ODD_NAME = "odd.mrcs"
EXT_MRC_EVEN_NAME = "even.mrc"
EXT_MRC_ODD_NAME = "odd.mrc"
EVEN = 'even'
ODD = 'odd'
MRCS_EXT = 'mrcs'
MRC_EXT = 'mrc'
EDF_EXT = 'edf'
XF_EXT = 'xf'
TLT_EXT = 'tlt'
RAWTLT_EXT = 'rawtlt'
DEFOCUS_EXT = 'defocus'
FID_EXT = 'fid'
TXT_EXT = 'txt'
REC_EXT = 'rec'
XYZ_EXT = 'xyz'
MOD_EXT = 'mod'
SEED_EXT = 'seed'
PREXF_EXT = 'prexf'
PREXG_EXT = 'prexg'
SFID_EXT = 'sfid'
RESID_EXT = 'resid'
