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
from typing import List, Set

import numpy as np

from pyworkflow.utils import cyanStr
from tomo.objects import TiltSeries, TiltImage

logger = logging.getLogger(__name__)


def genXfFile(ts: TiltSeries,
              outXfName: str,
              presentAcqOrders: Set[int] = None) -> None:
    """ This method takes a tilt series and the output transformation file path
    and creates an IMOD-based transform file in the location indicated. The transformation matrix
    of a tilt-image is only added if its acquisition order is contained in a set composed of the
    acquisition orders present in both the given tilt-series and the CTFTomoSeries. If presentAcqOrders
    is not None, it is considered before the attribute onlyEnabled, as presentAcqOrders may also have been
    generated considering the enabled elements of the intersection.
    """
    logger.info(cyanStr(f"tsId = {ts.getTsId()} -> Generating the transformation xf file {outXfName}..."))
    if presentAcqOrders:
        presentAcqOrders = ts.getTsPresentAcqOrders()
        tsMatrixList = [_tiMatrixToXfFormat(ti) for ti in ts if ti.getAcquisitionOrder() in presentAcqOrders]
    else:
        tsMatrixList = [_tiMatrixToXfFormat(ti) for ti in ts]

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


def fiducialModel2List(fiducialFilePath: str) -> list:
    """ This method takes an IMOD-based fiducial model file path
    and returns a list containing the coordinates of each
    fiducial for each tilt-image belonging to the tilt-series. """

    fiducialList = []

    with open(fiducialFilePath) as f:
        fiducialText = f.read().splitlines()

        for line in fiducialText:
            # Fix IMOD bug: columns merge in coordinates exceed 3 digits
            vector = line.replace('-', ' -').split()
            vector = [round(float(i)) for i in vector]
            fiducialList.append(vector)

    return fiducialList


def fidResidualModel2List(fiducialFilePath: str) -> list:
    """ This method takes an IMOD-based fiducial residual model
    file path and returns a list containing the coordinates
    and residual values of each fiducial for each tilt-image
    belonging to the tilt-series. Since IMOD establishes a
    float value for each coordinate they are parsed to int. """

    fiducialResidList = []

    with open(fiducialFilePath) as f:
        fiducialText = f.read().splitlines()

        for line in fiducialText[1:]:
            # Fix IMOD bug: columns merge in coordinates exceed 3 digits
            vector = line.replace('-', ' -').split()
            fiducialResidList.append([round(float(vector[0])),
                                      round(float(vector[1])),
                                      int(vector[2]),
                                      float(vector[3]),
                                      float(vector[4])])

    return fiducialResidList