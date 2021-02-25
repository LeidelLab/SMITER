======
SMITER
======


.. image:: https://img.shields.io/pypi/v/smiter.svg
        :target: https://pypi.python.org/pypi/smiter

.. image:: https://travis-ci.com/MKoesters/SMITER.svg?token=7Uh1o2G3gUZXxMB2xrqp&branch=dev
        :target: https://travis-ci.com/MKoesters/smiter

.. image:: https://readthedocs.org/projects/smiter/badge/?version=latest
        :target: https://smiter.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


.. image:: https://pyup.io/repos/github/MKoesters/smiter/shield.svg
     :target: https://pyup.io/repos/github/MKoesters/smiter/
     :alt: Updates


Summary
-------

Pytohn library to create synthetic mzMLs file based on chemical formulas. All molecules can be simulated du to abstraction to chemical formulas.

Abstract
-------

SMITER (Synthetic mzML writer) is a python-based command-line tool designed to simulate LC-MS/MS runs. It enables the simulation of any biomolecule since all calculations are based on the chemical formulas. As SMITER features a modular design, noise and fragmentation models can easily be implemented or adapted. By default, SMITER uses an established noise model and offers several methods for peptide fragmentation or two models for nucleoside fragmentation. Due to the rich python ecosystem, other modules, e.g. for retention time prediction, can easily be implemented for the tailored simulation of any molecule of choice. This allows for the facile creation of defined gold-standard-LC-MS/MS datasets for any type of experiment. Such gold standards, where the ground truth is known, are required in computational mass spectrometry to test new algorithms and to improve parameters for existing ones. Similarly, gold-standard datasets can be used to evaluate analytical hurdles e.g. by predicting co-elution and co-fragmentation of molecules. As these challenges hinder the detection or quantification of co-eluents, a comprehensive simulation can identify and thus prevent such difficulties before performing actual MS experiments. SMITER allows to create such datasets easily, fast and efficiently



* Free software: MIT license
* Documentation: https://smiter.readthedocs.io.


Features
--------

* TODO


Download and Installation
-------------------------

SMITER requires `Python`_ 3.6 or higher.


There are two recommended ways for installing SMITER:

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
~~~~~~~~~~~~~~~~~~~~~~~~

Clone the GitHub repo `GitHub`_::

   user@localhost:~$ git clone https://github.com/LeidelLal/SMITER.git


.. _GitHub:
   https://github.com/LeidelLab/SMITER


Install the requirements and SMITER::

    user@localhost:~$ cd smiter
    user@localhost:~/smiter$ pip install -r requirements.txt
    user@localhost:~/smiter$ python setup.py install

.. note::

	We recommedn using an virtual environment when using SMITER



Copyrights
***********

Copyright 2020-2021 by authors and contributors


* Manuel Koesters
* Johannes Leufken
* Sebastian Leidel


Contact
*******

    | Prof. Dr. Sebastian Leidel
    | University of Bern
	| Department of Chemistry, Biochemistry and Pharmaceutical Sciences
	| Freiestrasse 3
	| 3012 Bern
	| Switzerland


Citation
********

Please do not forget to cite SMITER:


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
