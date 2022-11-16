===========
IMOD plugin
===========

This plugin provides wrappers for several programs of `IMOD <https://bio3d.colorado.edu/imod/>`_ software suite.

.. image:: https://img.shields.io/pypi/v/scipion-em-imod.svg
        :target: https://pypi.python.org/pypi/scipion-em-imod
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-imod.svg
        :target: https://pypi.python.org/pypi/scipion-em-imod
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-imod.svg
        :target: https://pypi.python.org/pypi/scipion-em-imod
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-imod?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-imod
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-imod
        :target: https://pypi.python.org/pypi/scipion-em-imod
        :alt: Downloads

Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

   .. code-block::

      scipion installp -p scipion-em-imod

b) Developer's version

   * download repository

   .. code-block::

      git clone -b devel https://github.com/scipion-em/scipion-em-imod.git

   * install

   .. code-block::

      scipion installp -p /path/to/scipion-em-imod --devel

IMOD binaries will be downloaded and installed automatically with the plugin, but you can also link an existing installation. Default installation path assumed is ``software/em/imod-4.11.20/IMOD``, if you want to change it, set *IMOD_HOME* in ``scipion.conf`` file to the folder where the IMOD is installed.

To check the installation, simply run one of the tests. A complete list of tests can be displayed by executing ``scipion test --show --grep imod``

Supported versions
------------------

4.11.7, 4.11.20

Protocols
---------

* apply transformation
* automatic CTF estimation
* CTF correction
* coarse prealignment
* dose filter
* etomo interactive
* exclude views
* fiducial alignment
* generate fiducial model
* gold bead picker 3D
* import tomo CTFs
* import transformation matrix
* manual CTF estimation
* tilt-series preprocess
* tomo preprocess
* tomo projection
* tomo reconstruction
* X-rays eraser

References
----------

1. James R. Kremer, David N. Mastronarde, J.Richard McIntosh. Computer Visualization of Three-Dimensional Image Data Using IMOD. Journal of Structural Biology, Volume 116, Issue 1, 1996, Pages 71-76. https://doi.org/10.1006/jsbi.1996.0013
2. David N. Mastronarde, Susannah R. Held. Automated tilt series alignment and tomographic reconstruction in IMOD. Journal of Structural Biology, Volume 197, Issue 2, 2017, Pages 102-113, ISSN 1047-8477. https://doi.org/10.1016/j.jsb.2016.07.011
