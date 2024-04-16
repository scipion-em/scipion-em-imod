# *****************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
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

import pyworkflow.protocol.params as params


class CommonIMODforms:
    def trimimgForm(self, form, pxTrimCondition='False', correlationCondition='True', levelType=params.LEVEL_ADVANCED):
        '''
        Generally, this form will be integrated in a groupForm, the group form argument is form. A set of flags
        control what elements are shown
        '''
        form.addParam('pxTrim',
                          params.NumericListParam,
                          condition=pxTrimCondition,
                          label='Pixels to trim (x y without coma separator)',
                          default="0 0",
                          help='Pixels to trim off each side in X and Y.\n'
                               'Some trimming should be used for patch tracking',
                          expertLevel=levelType)

        xtrimming = form.addLine('Pixels to do correlation along X-axis',
                                     expertLevel=levelType,
                                     condition=correlationCondition,
                                     help="Starting and ending X coordinates of a region to correlate, "
                                          "based on the position of the region at zero tilt.")

        xtrimming.addParam('xmin',
                           params.IntParam,
                           label='X axis min (left)',
                           allowsNull=True,
                           expertLevel=levelType)

        xtrimming.addParam('xmax',
                           params.IntParam,
                           label='X axis max (right)',
                           allowsNull=True,
                           expertLevel=levelType)

        ytrimming = form.addLine('Pixels to do correlation along Y-axis',
                                     expertLevel=levelType,
                                     condition=correlationCondition,
                                     help="Starting and ending Y coordinates of a region to correlate, "
                                          "based on the position of the region at zero tilt.")

        ytrimming.addParam('ymin',
                           params.IntParam,
                           label='Y axis min (top)',
                           allowsNull=True,
                           expertLevel=levelType)

        ytrimming.addParam('ymax',
                           params.IntParam,
                           label='Y axis max (botton)',
                           allowsNull=True,
                           expertLevel=levelType)


    def filteringParametersForm(self, form, condition, levelType=params.LEVEL_NORMAL):
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
                                       ""
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
                       params.FloatParam,
                       label='Filter radius 1',
                       default='0.0',
                       expertLevel=levelType)

        line1.addParam('filterSigma1',
                       params.FloatParam,
                       label='Filter sigma 1',
                       default='0.03',
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
                       params.FloatParam,
                       label='Filter radius 2',
                       default='0.25',
                       expertLevel=levelType)

        line2.addParam('filterSigma2',
                       params.FloatParam,
                       label='Filter sigma 2',
                       default='0.05',
                       expertLevel=levelType)
