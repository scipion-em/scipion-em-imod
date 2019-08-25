# coding: latin-1
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
"""

@article{Kremer1996,
title = "Computer Visualization of Three-Dimensional Image Data Using IMOD",
journal = "Journal of Structural Biology",
volume = "116",
number = "1",
pages = "71 - 76",
year = "1996",
issn = "1047-8477",
doi = "https://doi.org/10.1006/jsbi.1996.0013",
url = "http://www.sciencedirect.com/science/article/pii/S1047847796900131",
author = "James R. Kremer and David N. Mastronarde and J.Richard McIntosh",
abstract = "We have developed a computer software package, IMOD, as a tool for analyzing and viewing three-dimensional biological image data. IMOD is useful for studying and modeling data from tomographic, serial section, and optical section reconstructions. The software allows image data to be visualized by several different methods. Models of the image data can be visualized by volume or contour surface rendering and can yield quantitative information."
}

@article{Mastronarde2017,
title = "Automated tilt series alignment and tomographic reconstruction in IMOD",
journal = "Journal of Structural Biology",
volume = "197",
number = "2",
pages = "102 - 113",
year = "2017",
note = "Electron Tomography",
issn = "1047-8477",
doi = "https://doi.org/10.1016/j.jsb.2016.07.011",
url = "http://www.sciencedirect.com/science/article/pii/S1047847716301526",
author = "David N. Mastronarde and Susannah R. Held",
keywords = "Electron tomography, Tilt series alignment, Tomographic reconstruction, Automated processing",
abstract = "Automated tomographic reconstruction is now possible in the IMOD software package, including the merging of tomograms taken around two orthogonal axes. Several developments enable the production of high-quality tomograms. When using fiducial markers for alignment, the markers to be tracked through the series are chosen automatically; if there is an excess of markers available, a well-distributed subset is selected that is most likely to track well. Marker positions are refined by applying an edge-enhancing Sobel filter, which results in a 20% improvement in alignment error for plastic-embedded samples and 10% for frozen-hydrated samples. Robust fitting, in which outlying points are given less or no weight in computing the fitting error, is used to obtain an alignment solution, so that aberrant points from the automated tracking can have little effect on the alignment. When merging two dual-axis tomograms, the alignment between them is refined from correlations between local patches; a measure of structure was developed so that patches with insufficient structure to give accurate correlations can now be excluded automatically. We have also developed a script for running all steps in the reconstruction process with a flexible mechanism for setting parameters, and we have added a user interface for batch processing of tilt series to the Etomo program in IMOD. Batch processing is fully compatible with interactive processing and can increase efficiency even when the automation is not fully successful, because users can focus their effort on the steps that require manual intervention."
}

"""
