======
SMITER
======


.. image:: https://img.shields.io/pypi/v/smiter.svg
        :target: https://pypi.python.org/pypi/smiter

.. image:: https://api.travis-ci.com/LeidelLab/SMITER.svg?branch=dev
        :target: https://travis-ci.com/LeidelLab/smiter

.. image:: https://readthedocs.org/projects/smiter/badge/?version=latest
        :target: https://smiter.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/LeidelLab/smiter/shield.svg
     :target: https://pyup.io/account/repos/github/LeidelLab/SMITER/
     :alt: Updates


Summary
-------

Python library to create synthetic mzMLs file based on chemical formulas. All molecules can be simulated du to abstraction to chemical formulas.

Abstract
-------

SMITER (Synthetic mzML writer) is a python-based command-line tool designed to simulate LC-MS/MS runs. It enables the simulation of any biomolecule since all calculations are based on the chemical formulas. As SMITER features a modular design, noise and fragmentation models can easily be implemented or adapted. By default, SMITER uses an established noise model and offers several methods for peptide fragmentation or two models for nucleoside fragmentation. Due to the rich python ecosystem, other modules, e.g. for retention time prediction, can easily be implemented for the tailored simulation of any molecule of choice. This allows for the facile creation of defined gold-standard-LC-MS/MS datasets for any type of experiment. Such gold standards, where the ground truth is known, are required in computational mass spectrometry to test new algorithms and to improve parameters for existing ones. Similarly, gold-standard datasets can be used to evaluate analytical hurdles e.g. by predicting co-elution and co-fragmentation of molecules. As these challenges hinder the detection or quantification of co-eluents, a comprehensive simulation can identify and thus prevent such difficulties before performing actual MS experiments. SMITER allows to create such datasets easily, fast and efficiently





Features
--------

* simulate mass spectrometry data for any biomolecule
* usage of highly-accurate isotopic patterns enabled by `pyQms`_
* feature scaling by gauss-, gamma- and exponentially-modified gaussian distributions
* m/z-and intensity noise injection ( uniform noise or a noise model that combines general noise with intensity-specific noise)
* MS2 fragmentation for peptides and modified nucleosides.
* Free software: MIT license
* Documentation: https://smiter.readthedocs.io.

.. _pyQms:
	https://github.com/pyQms/pyqms

Download and Installation
-------------------------

SMITER requires `Python`_ 3.7 or higher.


There are two recommended ways for installing SMITER

* Installation via pip
* Installation from the source (GitHub)

.. _Python:
   https://www.python.org/downloads/

.. _install_pip:

Installation via pip
--------------------

Execute the following command from your command line::

    user@localhost:~$ pip install smiter


Installation from source
------------------------

Clone the GitHub repo `GitHub`_::

   user@localhost:~$ git clone https://github.com/LeidelLab/SMITER.git


.. _GitHub:
   https://github.com/LeidelLab/SMITER


Install the requirements and SMITER::

    user@localhost:~$ cd smiter
    user@localhost:~/smiter$ pip install -r requirements.txt
    user@localhost:~/smiter$ python setup.py install


.. note::

	We recommend using an virtual environment when using SMITER



Testing
-------

To test the package and correct installation::

    user@localhost:~/smiter$ tox

Copyrights
----------

Copyright 2020-2021 by authors and contributors


* Manuel Kösters
* Johannes Leufken
* Sebastian Leidel


Contact
-------

 Prof. Dr. Sebastian Leidel
 University of Bern
 Department of Chemistry, Biochemistry and Pharmaceutical Sciences
 Freiestrasse 3
 3012 Bern
 Switzerland


Citation
--------

Please do not forget to cite SMITER:

<ref>


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
