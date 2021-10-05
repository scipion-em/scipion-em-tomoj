# **************************************************************************
# *
# * Authors:     Antoine Cossa (antoine.cossa@universite-paris-saclay.fr) [1]
# *
# * [1] Universite Paris-Saclay, Orsay, France
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
This module contains utils functions for TomoJ protocols
"""

import numpy as np


# def formatTransformFile(ts, transformFilePath):
#     """This method takes a tilt series and the output transformation file path
#     and creates an IMOD-based transform file in the location indicated"""
#     tsMatrixTransformList = []
#     for ti in ts:
#         transform = ti.getTransform().getMatrix().flatten()
#         transformIMOD = [transform[0],
#                          transform[1],
#                          transform[3],
#                          transform[4],
#                          transform[2],
#                          transform[5]]
#         tsMatrixTransformList.append(transformIMOD)
#     with open(transformFilePath, 'w') as f:
#         csvW = csv.writer(f, delimiter='\t')
#         csvW.writerows(tsMatrixTransformList)


def formatTransformationMatrix(matrixFile):
    """This method takes the TomoJ-based transformation matrix file path
    and returns a 3D matrix containing the transformation matrices for each
    tilt-image belonging to the tilt-series.
    TomoJ transformation files are formatted as follows:
    [tilt angle][tilt axis][0,0][1,0][0,1][1,1][0,2][1,2]"""
    with open(matrixFile, "r") as matrix:
        lines = matrix.readlines()
    numberLines = len(lines)
    frameMatrix = np.empty([3, 3, numberLines])
    i = 0
    # for line in lines:
    #     print(line)
    #     values = line.split()
    #     frameMatrix[0, 0, i] = float(values[2])
    #     frameMatrix[1, 0, i] = float(values[3])
    #     frameMatrix[0, 1, i] = float(values[4])
    #     frameMatrix[1, 1, i] = float(values[5])
    #     frameMatrix[0, 2, i] = float(values[6])
    #     frameMatrix[1, 2, i] = float(values[7])
    #     frameMatrix[2, 0, i] = 0.0
    #     frameMatrix[2, 1, i] = 0.0
    #     frameMatrix[2, 2, i] = 1.0
    #     i += 1

    # DEBUG version
    for line in lines:
        if i != numberLines-1:
            values = line.split()
            frameMatrix[0, 0, i] = float(values[2])
            frameMatrix[1, 0, i] = float(values[3])
            frameMatrix[0, 1, i] = float(values[4])
            frameMatrix[1, 1, i] = float(values[5])
            frameMatrix[0, 2, i] = float(values[6])
            frameMatrix[1, 2, i] = float(values[7])
            frameMatrix[2, 0, i] = 0.0
            frameMatrix[2, 1, i] = 0.0
            frameMatrix[2, 2, i] = 1.0
            i += 1
    return frameMatrix