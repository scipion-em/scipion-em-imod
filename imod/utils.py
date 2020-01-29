# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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
"""
This module contains utils functions for IMOD protocols
"""

import csv

def formatTransformFile(ts, transformFilePath):
    tsMatrixTransformList = []
    for ti in ts:
        transform = ti.getTransform().getMatrix().flatten()
        transformIMOD = [transform[0],
                         transform[1],
                         transform[3],
                         transform[4],
                         transform[2],
                         transform[5]]
        tsMatrixTransformList.append(transformIMOD)
    with open(transformFilePath, 'w') as f:
        csvW = csv.writer(f, delimiter='\t')
        csvW.writerows(tsMatrixTransformList)
