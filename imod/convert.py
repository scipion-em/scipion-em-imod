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
import csv
import logging
from typing import List

import numpy as np

from pyworkflow.utils import cyanStr
from tomo.objects import TiltSeries, TiltImage

logger = logging.getLogger(__name__)


def genTltFile(ts: TiltSeries,
               outTltFileName: str,
               ignoreExcludedViews: bool = False) -> None:
    logger.info(cyanStr(f"tsId = {ts.getTsId()} -> Generating the angles tlt file {outTltFileName}..."))
    ts.generateTltFile(outTltFileName,
                       excludeViews=not ignoreExcludedViews)


def genXfFile(ts: TiltSeries,
              outXfName: str,
              ignoreExcludedViews: bool = False) -> None:
    """ This method takes a tilt series and the output transformation file path
    and creates an IMOD-based transform file in the location indicated. The transformation matrix
    of a tilt-image is only added if its acquisition order is contained in a set composed of the
    acquisition orders present in both the given tilt-series and the CTFTomoSeries. If presentAcqOrders
    is not None, it is considered before the attribute onlyEnabled, as presentAcqOrders may also have been
    generated considering the enabled elements of the intersection.
    """
    logger.info(cyanStr(f"tsId = {ts.getTsId()} -> Generating the transformation xf file {outXfName}..."))
    if ignoreExcludedViews:
        tsMatrixList = [_tiMatrixToXfFormat(ti) for ti in ts]
    else:
        presentAcqOrders = ts.getTsPresentAcqOrders()
        tsMatrixList = [_tiMatrixToXfFormat(ti) for ti in ts if ti.getAcquisitionOrder() in presentAcqOrders]

    with open(outXfName, 'w') as f:
        csvW = csv.writer(f, delimiter='\t')
        csvW.writerows(tsMatrixList)

def _tiMatrixToXfFormat(tiltImage: TiltImage) -> List[str]:
    transform = tiltImage.getTransform().getMatrix().flatten()
    transformIMOD = ['%.7f' % transform[0],
                     '%.7f' % transform[1],
                     '%.7f' % transform[3],
                     '%.7f' % transform[4],
                     "{:>6}".format('%.3g' % transform[2]),
                     "{:>6}".format('%.3g' % transform[5])]
    return transformIMOD


def readXfFile(xfFile) -> np.ndarray:
    """ This method takes an IMOD-based transformation matrix file (.xf) path and
    returns a 3D matrix containing the transformation matrices for
    each tilt-image belonging to the tilt-series. """

    matrix = np.loadtxt(xfFile, dtype=float, comments='#')
    numberLines = matrix.shape[0]
    frameMatrix = np.empty([3, 3, numberLines])

    for row in range(numberLines):
        frameMatrix[0, 0, row] = matrix[row][0]
        frameMatrix[1, 0, row] = matrix[row][2]
        frameMatrix[0, 1, row] = matrix[row][1]
        frameMatrix[1, 1, row] = matrix[row][3]
        frameMatrix[0, 2, row] = matrix[row][4]
        frameMatrix[1, 2, row] = matrix[row][5]
        frameMatrix[2, 0, row] = 0.0
        frameMatrix[2, 1, row] = 0.0
        frameMatrix[2, 2, row] = 1.0

    return frameMatrix