# **************************************************************************
# *
# * Authors:     Antoine Cossa (antoine.cossa@universite-paris-saclay.fr) [1]
#                Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [2]
# *
# * [1] Universite Paris-Saclay, Orsay, France
# * [2] Centro Nacional de Biotecnologia, CSIC, Spain
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
import pwem

from .constants import TOMOJ_HOME, DEFAULT_VERSION

__version__ = '3.0.8'
_logo = "tomoj_icon.png"
_references = ['MessaoudiI2007', 'Sorzano2009']


class Plugin(pwem.Plugin):
    _homeVar = TOMOJ_HOME
    _supportedVersions = DEFAULT_VERSION

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(TOMOJ_HOME, cls._getTOMOJFolder(DEFAULT_VERSION))

    @classmethod
    def _getEMFolder(cls, version, *paths):
        return os.path.join("tomoj-%s" % version, *paths)

    @classmethod
    def _getTOMOJFolder(cls, version, *paths):
        return os.path.join(cls._getEMFolder(version), *paths)

    @classmethod
    def getEnviron(cls):
        env = pwem.pwutils.Environ(os.environ)
        if 'TOMOJ_DIR' in env:
            del env['TOMOJ_DIR']
        if 'TOMOJ_PATH' in env:
            del env['TOMOJ_PATH']
        return env

    @classmethod
    def getDependencies(cls):
        neededPrograms = ['java', 'python']

        return neededPrograms

    @classmethod
    def defineBinaries(cls, env):
        TOMOJ_INSTALLED = 'tomoj_%s_installed' % DEFAULT_VERSION

        # Download TomoJ .jar
        # installationCmd = 'wget --continue --https-only ' \
        #                   'https://sourceforge.net/projects/tomoj/files/' \
        #                   'TomoJ_Applications-%s-jar-with-dependencies.jar && ' \
        #                   % DEFAULT_VERSION
        # Temporary solution to get TomoJ's latest version
        installationCmd = 'wget --continue ' \
                          'http://xfer.curie.fr/get/KaAoWTgodnt/' \
                          'TomoJ_Applications-%s-jar-with-dependencies.jar && '\
                          % DEFAULT_VERSION

        # Install dependencies
        # Fractional_Splines_Wavelets.jar
        installationCmd += 'wget --continue ' \
                           'http://bigwww.epfl.ch/demo/fractsplines/java/' \
                           'Fractional_Splines_Wavelets.jar && '

        # Create installation finished flag file
        installationCmd += 'touch %s' % TOMOJ_INSTALLED

        env.addPackage('tomoj',
                       version=DEFAULT_VERSION,
                       tar='void.tgz',
                       createBuildDir=True,
                       buildDir=pwem.Config.EM_ROOT,
                       neededProgs=cls.getDependencies(),
                       commands=[(installationCmd, TOMOJ_INSTALLED)],
                       default=True)

    @classmethod
    def runTomoJ(cls, protocol, args, cwd=None):
        """ Run TomoJ command from a given protocol. """

        # Get the command
        cmd = "java -cp " + os.path.join(pwem.Config.EM_ROOT,
                                         cls._getTOMOJFolder(DEFAULT_VERSION))\
              + "/TomoJ_Applications-%s-jar-with-dependencies.jar" \
              % DEFAULT_VERSION + " fr.curie.tomoj.TomoJ "

        # Run the protocol with that command
        protocol.runJob(cmd, args, env=cls.getEnviron(), cwd=cwd)

