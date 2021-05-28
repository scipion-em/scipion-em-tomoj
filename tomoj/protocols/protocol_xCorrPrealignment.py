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

import os
import numpy as np
import imod.utils as utils
import pwem.objects as data
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase


class ProtTomojXcorrPrealignment(EMProtocol, ProtTomoBase):
    """
    Tilt-series' cross correlation alignment based on the TomoJ procedure.
    More info:
        http://u759.sfbiophys.org/software/update/20140207/Manual_TomoJ_2.24.pdf DEAD LINK
    """

    _label = 'xcorr prealignment'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-Series.')

        form.addParam('integerTranslation', params.BooleanParam,
                      default=False,
                      label='Compute integer translation',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Compute integer translation via cross correlation '
                           'for the tilt-series.')

        form.addParam('downsampling', params.IntParam,
                      default=1,
                      label='Downsampling',
                      help='Reduce the size of images for computation. The '
                           'value (integer) corresponds to the factor of '
                           'reduction, usually 2,4,8...')

        form.addParam('roi', params.BooleanParam,
                      default=False,
                      label='ROI centered cross-correlation',
                      important=True,
                      help='Take the central part of images of size rwidth '
                           'rheight (integers) to compute cross-correlation.')

        roi = form.addLine('', condition='roi',
                           help='Width and height of the centered ROI.')

        roi.addParam('rwidth', params.IntParam,
                     label='Width',
                     important=True,
                     help='Width of the centered ROI.')

        roi.addParam('rheight', params.IntParam,
                     label='Height',
                     important=True,
                     help='Height of the centered ROI.')

        bandpass = form.addLine('Bandpass filter',
                                help='Apply a bandpassfilter on images (after '
                                     'roi and downsampling if any). The 2 '
                                     'first values (double) correspond to the '
                                     'radius of the band in pixels [minimum '
                                     'maximum]. The third value corresponds '
                                     '(double) to the sinusoidal decrease '
                                     'radius (in pixels) to prevent artifacts.')

        bandpass.addParam('bandpassmin', params.FloatParam,
                          label='Min',
                          important=True)

        bandpass.addParam('bandpassmax', params.FloatParam,
                          label='Max',
                          important=True)
        bandpass.addParam('bandpassdecrease', params.IntParam,
                          label='Decrease',
                          important=True)

        form.addParam('variancefilter', params.IntParam,
                      default=1,
                      label='Variance filter radius',
                      important=True,
                      help='Apply a variance filter on images with the given '
                           'radius (integer). It results in contours images.')
        form.addParam('expandimage', params.FloatParam,
                      default=0.0,
                      label='Expand images',
                      important=True,
                      help='Expands the image to correct the stretching due to'
                           ' tilt. To do this correctly the tilt axis needs to'
                           ' be given as angle (double) from vertical axis.')

        form.addParam('multiscale', params.IntParam,
                      default=1,
                      label='Multiscale cross-correlation',
                      important=True,
                      help='Apply a multiscale approach with the given number '
                           'of level.')

        form.addParam('cumulativereference', params.BooleanParam,
                      default=False,
                      label='Cumulative reference',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='If true, the processing is not done between '
                           'consecutive images but using central image as '
                           'reference to which is added the newly aligned '
                           'images.')

        form.addParam('loop', params.BooleanParam,
                      default=False,
                      label='Loop until stabilization',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Refine alignment by doing the cross-correlation as'
                           ' many times as needed to stabilize.')

        form.addParam('computeAlignment', params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Generate interpolated tilt-series', important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Generate and save the interpolated tilt-series '
                           'applying the obtained transformation matrices.')

        group = form.addGroup('Interpolated tilt-series',
                              condition='computeAlignment==0')

        group.addParam('binning', params.FloatParam,
                       default=1.0,
                       label='Binning',
                       help='Binning to be applied to the interpolated '
                            'tilt-series. Must be a integer bigger than 1.')

        # form.addParam('rotationAngle',
        #               params.FloatParam,
        #               label='Tilt rotation angle (deg)',
        #               default='0.0',
        #               expertLevel=params.LEVEL_ADVANCED,
        #               help="Angle from the vertical to the tilt axis in raw images.")

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('computeXcorrStep', ts.getObjId())
            # if self.computeAlignment.get() == 0:
            #     self._insertFunctionStep('computeInterpolatedStackStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)
        outputTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)

        """Apply the transformation from the input tilt-series"""
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        ts.generateTltFile(angleFilePath)

    def computeXcorrStep(self, tsObjId):
        """Compute transformation matrix for each tilt series"""
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        # integerTranslation = ""
        # if self.integerTranslation: integerTranslation = 'integertranslation'
        #
        # cumulativereference = ""
        # if self.cumulativereference: cumulativereference = 'cumulativereference'
        #
        # loop = ""
        # if self.loop.get() == True: loop = 'loop'

        paramsXcorr = {
            'input': os.path.join(tmpPrefix, '%s.st' % tsId),
            'output': os.path.join(extraPrefix, '%s.prexf' % tsId),
            'tiltfile': os.path.join(tmpPrefix, '%s.rawtlt' % tsId),
            # 'integerTranslation': integerTranslation,
            'downsampling': self.downsampling.get(),
            'rwidth': self.rwidth.get(),
            'rheight': self.rheight.get(),
            'bandpassmin': self.bandpassmin.get(),
            'bandpassmax': self.bandpassmax.get(),
            'bandpassdecrease': self.bandpassdecrease.get(),
            'variancefilter': self.variancefilter.get(),
            'expandimage': self.expandimage.get(),
            'multiscale': self.multiscale.get(),
            # 'cumulativereference': cumulativereference,
            # 'loop': loop
            # 'RotationAngle': self.rotationAngle.get(),
            # 'FilterSigma1': 0.03,
            # 'FilterSigma2': 0.05,
            # 'FilterRadius2': 0.25
        }
        # argsXcorr = "-loadangles %(tiltfile)s " \
        #             "-xcorr " \
        #             "%(integerTranslation)s " \
        #             "downsampling %(downsampling)f " \
        #             "roi %(rwidth)d %(rheight)d " \
        #             "bandpassfilter %(bandpassmin)f %(bandpassmax)f %(bandpassdecrease)d" \
        #             "variancefilter %(variancefilter)d " \
        #             "expandimage %(expandimage)f " \
        #             "multiscale %(multiscale)d " \
        #             "%(cumulativereference)s " \
        #             "%(loop)s " \
        #             "/home/acossa/ScipionUserData/projects/TestImodReconstructionWorkflow/%(input)s "  # Input is the last argument

        argsXcorr = "-loadangles %(tiltfile)s " \
                    "-xcorr " \
                    "downsampling %(downsampling)d " \
                    "roi %(rwidth)d %(rheight)d " \
                    "bandpassfilter %(bandpassmin)f %(bandpassmax)f %(bandpassdecrease)d " \
                    "variancefilter %(variancefilter)d " \
                    "expandimage %(expandimage)f " \
                    "multiscale %(multiscale)d "
                    # "-output %(output)s " \
                    # "-RotationAngle %(RotationAngle)f " \
                    # "-FilterSigma1 %(FilterSigma1)f " \
                    # "-FilterSigma2 %(FilterSigma2)f " \
                    # "-FilterRadius2 %(FilterRadius2)f"
        if self.integerTranslation:
            argsXcorr += 'integertranslation '
        if self.cumulativereference:
            argsXcorr += 'cumulativereference '
        if self.loop:
            argsXcorr += 'loop '

        # Add input as last arg
        argsXcorr += "/home/acossa/ScipionUserData/projects/TestImodReconstructionWorkflow/%(input)s "
        print(argsXcorr)
        self.runJob('/home/acossa/ImageJ/jre/bin/java -Xmx28000m -cp /home/acossa/ImageJ/plugins/TomoJ_Applications-2.7-jar-with-dependencies.jar fr.curie.tomoj.TomoJ ', argsXcorr % paramsXcorr)




        # paramsXftoxg = {
        #     'input': os.path.join(extraPrefix, '%s.prexf' % tsId),
        #     #'output': os.path.join(extraPrefix, '%s.prexg' % tsId),
        # }
        # argsXftoxg = "-input %(input)s " \
        #              "-goutput %(goutput)s"
        # self.runJob('xftoxg', argsXftoxg % paramsXftoxg)

        """Generate output tilt series"""
        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()
        tsId = ts.getTsId()
        alignmentMatrix = utils.formatTransformationMatrix(self._getExtraPath('%s/%s.prexg' % (tsId, tsId)))
        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputSetOfTiltSeries.append(newTs)
        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(tiltImage.getLocation())
            transform = data.Transform()
            transform.setMatrix(alignmentMatrix[:, :, index])
            newTi.setTransform(transform)
            newTs.append(newTi)
        newTs.write()
        outputSetOfTiltSeries.update(newTs)
        outputSetOfTiltSeries.write()
        self._store()

    def computeInterpolatedStackStep(self, tsObjId):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()
        ts = self.inputSetOfTiltSeries.get()[tsObjId]

        tsId = ts.getTsId()
        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputInterpolatedSetOfTiltSeries.append(newTs)
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsAlignment = {
            'input': os.path.join(tmpPrefix, '%s.st' % tsId),
            'output': os.path.join(extraPrefix, '%s_preali.st' % tsId),
            'xform': os.path.join(extraPrefix, "%s.prexg" % tsId),
            'bin': int(self.binning.get()),
            'imagebinned': 1.0
        }
        argsAlignment = "-input %(input)s " \
                        "-output %(output)s " \
                        "-xform %(xform)s " \
                        "-bin %(bin)d " \
                        "-imagebinned %(imagebinned)s"
        self.runJob('newstack', argsAlignment % paramsAlignment)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1, (os.path.join(extraPrefix, '%s_preali.st' % tsId)))
            if self.binning > 1:
                newTi.setSamplingRate(tiltImage.getSamplingRate() * int(self.binning.get()))
            newTs.append(newTi)
        if self.binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))
        newTs.write()
        outputInterpolatedSetOfTiltSeries.update(newTs)  # update items and size info
        outputInterpolatedSetOfTiltSeries.write()
        self._store()

        """Debug code"""
        path.moveTree(self._getTmpPath(), self._getExtraPath())

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfTiltSeries(self):
        if not hasattr(self, "outputSetOfTiltSeries"):
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    def getOutputInterpolatedSetOfTiltSeries(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)
        return self.outputInterpolatedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfTiltSeries.getSize()))
        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           "Interpolated Tilt-Series: %d.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           % (self.outputSetOfTiltSeries.getSize()))
        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           "Also, interpolation has been completed for %d Tilt-series.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
