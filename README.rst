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

Current development
-------------------

This plugin is currently in **BETA** mode.


Installation
------------

You will need to use `3.0+ <https://scipion-em.github.io/docs/release-3.0.0/docs/scipion-modes/how-to-install.html>`_ version of Scipion to be able to run these protocols.

Protocols
---------

* **Apply transformation** : Compute the interpolated tilt-series from its transform matrix. More info: `newstack doc <https://bio3d.colorado.edu/imod/doc/man/newstack.html>`__
* **Automatic CTF estimation** :  CTF estimation of a set of input tilt-series using the IMOD procedure. More info: `ctfplotter doc <https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html>`_
* **CTF correction** : CTF correction of a set of input tilt-series using the IMOD procedure. More info: `ctfphaseflip doc <https://bio3d.colorado.edu/imod/doc/man/ctfphaseflip.html>`_
* **Dose filter** : Tilt-series dose filtering based on the IMOD procedure. More info: `mtffilter doc <https://bio3d.colorado.edu/imod/doc/man/mtffilter.html>`_
* **Etomo interactive** : Simple wrapper around etomo to manually reconstruct a Tomogram. More info:  `etomo tutorial <https://bio3d.colorado.edu/imod/doc/etomoTutorial.html>`_
* **Exclude views** : excludeviews - Reversibly remove views from a tilt series stack. If you use this protocol, make sure tis output tilt series is use for everything else  CTF estimation, per particle per tilt, tomogram reconstruction....More info:  `here <https://bio3d.colorado.edu/imod/doc/man/excludeviews.html>`_
* **Fiducial alignment** : Construction of a fiducial model and alignment of tilt-series based on the IMOD procedure. More info: `tiltalign doc <https://bio3d.colorado.edu/imod/doc/man/tiltalign.html>`_ , `model2point doc <https://bio3d.colorado.edu/imod/doc/man/model2point.html>`_, `imodtrans doc <https://bio3d.colorado.edu/imod/doc/man/imodtrans.html>`_, `newstack doc <https://bio3d.colorado.edu/imod/doc/man/newstack.html>`__, `ccderaser doc <https://bio3d.colorado.edu/imod/doc/man/ccderaser.html>`_
* **Generate fiducial model** : Construction of a fiducial model and alignment of tilt-series based on the IMOD procedure. More info: `autofidseed doc <https://bio3d.colorado.edu/imod/doc/man/autofidseed.html>`_, `beadtrack doc <https://bio3d.colorado.edu/imod/doc/man/beadtrack.html>`_, `model2point doc <https://bio3d.colorado.edu/imod/doc/man/model2point.html>`_
* **Gold bead picker 3D** : 3-dimensional gold bead picker using the IMOD procedure. More info: `findbeads3d doc <https://bio3d.colorado.edu/imod/doc/man/findbeads3d.html>`_
* **Import tomo CTFs** :  Protocol to import estimations of CTF series from tilt-series into Scipion.
* **Import transformation matrix** : Import the transformation matrices assigned to an input set of tilt-series
* **Manual CTF estimation** : CTF estimation of a set of input tilt-series using the IMOD procedure. Runs the protocol through the interactive GUI. The resulting defocus values MUST BE SAVED manually by the user. More info: `ctfplotter doc <https://bio3d.colorado.edu/imod/doc/man/ctfplotter.html>`_
* **Tilt-series preprocess** : Normalize input tilt-series and change its storing formatting. More info: `newstack doc <https://bio3d.colorado.edu/imod/doc/man/newstack.html>`__
* **Tomo preprocess** : Normalize input tomogram and change its storing formatting. More info: `newstack doc <https://bio3D.colorado.edu/imod/doc/newstack.html>`__, `binvol doc <https://bio3D.colorado.edu/imod/doc/binvol.html>`_
* **Tomo projection** : Re-project a tomogram given a geometric description (axis and angles). More info: `xyzproj doc <https://bio3d.colorado.edu/imod/doc/man/xyzproj.html>`_
* **Tomo reconstruction** : omogram reconstruction procedure based on the IMOD procedure. More info: `tilt doc <https://bio3d.colorado.edu/imod/doc/man/tilt.html>`_
* **Coarse prealignment** : Tilt-series cross correlation alignment based on the IMOD procedure. More info: `tiltxcorr doc <https://bio3d.colorado.edu/imod/doc/man/tiltxcorr.html>`_
* **X-rays eraser** : Erase X-rays from aligned tilt-series based on the IMOD procedure. More info: `ccderaser doc <https://bio3d.colorado.edu/imod/doc/man/ccderaser.html>`_

**Latest plugin versions**
==========================

If you want to check the latest version and release history go to `CHANGES <https://github.com/scipion-em/scipion-em-imod/imod/blob/master/CHANGES.txt>`_


**Installing the plugin**
=========================

In order to install the plugin follow these instructions:

.. code-block::

    scipion installp -p scipion-em-imod


or through the **plugin manager** by launching Scipion and following **Configuration** >> **Plugins**


**To install in development mode**

Clone or download the plugin repository

.. code-block::

    git clone https://github.com/scipion-em/scipion-em-imod.git

Install the plugin in developer mode.

.. code-block::

    scipion installp -p local/path/to/scipion-em-imod --devel


IMOD binaries will be downloaded and installed automatically with the plugin, but you can also link an existing installation. Default installation path assumed is ``software/em/imod-4.11.24/IMOD``, if you want to change it, set *IMOD_HOME* in ``scipion.conf`` file to the folder where the IMOD is installed.

To check the installation, simply run one of the tests. A complete list of tests can be displayed by executing ``scipion test --show --grep imod``

Supported versions
------------------

4.11.7, 4.11.20, 4.11.24


References
----------

1. James R. Kremer, David N. Mastronarde, J.Richard McIntosh. Computer Visualization of Three-Dimensional Image Data Using IMOD. Journal of Structural Biology, Volume 116, Issue 1, 1996, Pages 71-76. https://doi.org/10.1006/jsbi.1996.0013
2. David N. Mastronarde, Susannah R. Held. Automated tilt series alignment and tomographic reconstruction in IMOD. Journal of Structural Biology, Volume 197, Issue 2, 2017, Pages 102-113, ISSN 1047-8477. https://doi.org/10.1016/j.jsb.2016.07.011


Buildbot status
---------------

Status devel version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/imod_devel.svg


Status production version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/imod_prod.svg

