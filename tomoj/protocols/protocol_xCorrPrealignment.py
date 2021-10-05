# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
# *              Antoine Cossa (antoine.cossa@universite-paris-saclay.fr) [2]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# * [2] Universite Paris-Saclay, Orsay, France
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
import utils
import pwem.objects as data
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow import BETA
from pyworkflow.object import Set
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from tomoj import Plugin


class ProtTomojXcorrPrealignment(EMProtocol, ProtTomoBase):
    """
    Tilt-series' cross correlation alignment based on the TomoJ procedure.
    """

    _label = 'xcorr prealignment'
    _devStatus = BETA

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
                      default=True,
                      label='Compute integer translation',
                      important=True,
                      help='Compute integer translation via cross correlation '
                           'for the tilt-series.')

        form.addParam('downsampling', params.IntParam,
                      default=1,
                      label='Downsampling',
                      help='Reduce the size of images for computation. The '
                           'value (integer) corresponds to the factor of '
                           'reduction, usually 2,4,8...',
                      expertLevel=params.LEVEL_ADVANCED)

        roi = form.addLine('ROI centered cross-correlation',
                           help='Take the central part of images of size Width'
                                ' Height (integers) to compute cross-'
                                'correlation.',
                           expertLevel=params.LEVEL_ADVANCED)

        roi.addParam('roi', params.BooleanParam,
                     default=False,
                     label='',
                     important=True)

        roi.addParam('rwidth', params.IntParam,
                     label='Width',
                     important=True,
                     condition='roi',
                     help='Width of the centered ROI.')

        roi.addParam('rheight', params.IntParam,
                     label='Height',
                     important=True,
                     condition='roi',
                     help='Height of the centered ROI.')

        bandpass = form.addLine('Bandpass filter',
                                help='Apply a bandpassfilter on images (after '
                                     'roi and downsampling if any). The 2 '
                                     'first values (double) correspond to the '
                                     'radius of the band in pixels [minimum '
                                     'maximum]. The third value corresponds '
                                     '(double) to the sinusoidal decrease '
                                     'radius (in pixels) to prevent artifacts.',
                                expertLevel=params.LEVEL_ADVANCED)

        bandpass.addParam('bandpass', params.BooleanParam,
                          default=False,
                          important=True)

        bandpass.addParam('bandpassmin', params.FloatParam,
                          label='Min',
                          important=True,
                          condition='bandpass')

        bandpass.addParam('bandpassmax', params.FloatParam,
                          label='Max',
                          important=True,
                          condition='bandpass')

        bandpass.addParam('bandpassdecrease', params.IntParam,
                          label='Decrease',
                          important=True,
                          condition='bandpass')

        form.addParam('variancefilter', params.IntParam,
                      default=1,
                      label='Variance filter radius',
                      important=True,
                      help='Apply a variance filter on images with the given '
                           'radius (integer). It results in contours images.',
                      expertLevel=params.LEVEL_ADVANCED)

        expand = form.addLine('Expand images',
                              help='Expands the image to correct the '
                                   'stretching due to tilt. To do this '
                                   'correctly the tilt axis needs to be given '
                                   'as angle (double) from vertical axis.',
                              expertLevel=params.LEVEL_ADVANCED)

        expand.addParam('expand', params.BooleanParam,
                        default=False,
                        important=True)

        expand.addParam('expandimage', params.FloatParam,
                        default=0.0,
                        label='Tilt-axis',
                        important=True,
                        condition='expand')

        form.addParam('multiscale', params.IntParam,
                      default=2,
                      label='Multiscale cross-correlation',
                      important=True,
                      help='Apply a multiscale approach with the given number '
                           'of level. (1 = no multiscale)',
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('cumulativereference', params.BooleanParam,
                      default=False,
                      label='Cumulative reference',
                      important=True,
                      help='If true, the processing is not done between '
                           'consecutive images but using central image as '
                           'reference to which is added the newly aligned '
                           'images.',
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('loop', params.BooleanParam,
                      default=True,
                      label='Loop until stabilization',
                      important=True,
                      help='Refine alignment by doing the cross-correlation as'
                           ' many times as needed to stabilize.')

        form.addParam('computeAlignment', params.BooleanParam,
                      default=False,
                      label='Generate interpolated tilt-series',
                      important=True,
                      help='Generate and save the interpolated tilt-series '
                           'applying the obtained transformation matrices.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('computeXcorrStep', ts.getObjId())
        self._insertFunctionStep('closeOutputSetsStep')

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        """Apply the transformation from the input tilt-series"""
        outputTsFileName = os.path.join(tmpPrefix,
                                        ts.getFirstItem().parseFileName())
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(
            tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt"))
        ts.generateTltFile(angleFilePath)

    def computeXcorrStep(self, tsObjId):
        """Compute transformation matrix for each tilt-series"""
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsXcorr = {
            # 'input': os.path.join(tmpPrefix, '%s.st' % tsId),
            'input': os.path.join(tmpPrefix, ts.getFirstItem().parseFileName()),
            # 'output': os.path.join(extraPrefix, '%s.prexf' % tsId),
            # 'tiltfile': os.path.join(tmpPrefix, '%s.rawtlt' % tsId),
            'tiltfile': os.path.join(
                tmpPrefix, ts.getFirstItem().parseFileName(extension=".tlt")),
            'downsampling': self.downsampling.get(),
            'variancefilter': self.variancefilter.get(),
            'multiscale': self.multiscale.get(),
        }

        argsXcorr = "-loadangles %(tiltfile)s " \
                    "-xcorr " \
                    "downsampling %(downsampling)d " \
                    "variancefilter %(variancefilter)d " \
                    "multiscale %(multiscale)d "

        if self.integerTranslation:
            argsXcorr += 'integertranslation '
        if self.cumulativereference:
            argsXcorr += 'cumulativereference '
        if self.loop:
            argsXcorr += 'loop '
        if self.roi:
            argsXcorr += 'roi %d %d ' % (self.rwidth.get(), self.rheight.get())
        if self.bandpass:
            argsXcorr += 'bandpassfilter %f %f %d ' % (self.bandpassmin.get(),
                                                       self.bandpassmax.get(),
                                                       self.bandpassdecrease.get())
        if self.expand:
            argsXcorr += 'expandimage %f ' % self.expandimage.get()
        if self.computeAlignment:
            argsXcorr += '-savealignedimages '

        # Add input as last arg
        argsXcorr += "%(input)s "

        # Run TomoJ
        Plugin.runTomoJ(self, argsXcorr % paramsXcorr)

        """Debug code"""
        # path.moveTree(self._getTmpPath(), self._getExtraPath())
        """Move files to extra path"""
        if self.computeAlignment:
            path.moveFile(os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(
                              suffix="_ali", extension=".mrc")),
                          os.path.join(extraPrefix, ts.getFirstItem().parseFileName(
                              extension=".mrc")))

        path.moveFile(os.path.join(tmpPrefix, ts.getFirstItem().parseFileName(
                          extension="_xcorr.txt")),
                      os.path.join(extraPrefix, ts.getFirstItem().parseFileName(
                          extension="_xcorr.txt")))

        """Generate output tilt series"""
        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()
        alignmentMatrix = utils.formatTransformationMatrix(
            os.path.join(extraPrefix, ts.getFirstItem().parseFileName(
                extension="_xcorr.txt")))
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

        newTs.write(properties=False)

        outputSetOfTiltSeries.update(newTs)
        outputSetOfTiltSeries.write()

        self._store()

    def closeOutputSetsStep(self):
        self.getOutputSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfTiltSeries(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
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
        if not self.computeAlignment:
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfTiltSeries.getSize()))
        elif self.computeAlignment:
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           "Interpolated Tilt-Series: %d.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputSetOfTiltSeries.getSize(),
                              self.outputSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if not self.computeAlignment:
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the TomoJ procedure.\n"
                           % (self.outputSetOfTiltSeries.getSize()))
        elif self.computeAlignment:
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the TomoJ procedure.\n"
                           "Also, interpolation has been completed for %d Tilt-series.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
