# GDMATE - GeoDynamic Modeling Analysis Toolkit and Education
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gdmate/gdmate/HEAD)

## About

GDMATE is a Python library for generating and analyzing geodynamic models, along with related education.

Documentation: [http://gdmate.readthedocs.io](http://gdmate.readthedocs.io)

Source code: [https://github.com/gdmate/gdmate](https://github.com/gdmate/gdmate)

Authors (as of 2022)
* Dylan Vasey
* John Naliboff

## Requirements

* Python 3.7+
* Python modules:
  NumPy, SciPy, Matplotlib, Pyvista

## Installation

GDMATE is still in the earliest stages of development and is not yet available from typical hosting platforms (PyPI, Anaconda, etc.). For the moment, there are 2 ways to use GDMATE.

1. Click on the Binder badge at the top of this README to launch a Python environment with GDMATE installed in your web browser.

2. Clone this repository, then install GDMATE into a Python environment (using `virutalenv` or `conda`) using `pip`

```
git clone https://github.com/gdmate/gdmate.git
cd gdmate
pip install .
```
This should install GDMATE and its dependencies into your current environment. If you aren't familiar with managing Python virtual environments, [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-envs) is a good place to start.

## Getting started: The GDMATE Notebooks

For the moment, the main features of GDMATE are illustrated in the suite of Jupyter Notebooks housed in the `notebooks` directory. You can run and modifty these in the Binder environment [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gdmate/gdmate/HEAD) or in your local environment with GDMATE installed. Additional notebooks and example scripts will be added as development proceeds.

## About scripting in Python

GDMATE has the advantage of being adaptable and extensible in easy scripts.
As GDMATE is a toolkit, a graphical user interface would be impractical.
Nevertheless, we hope that we have succeeded in making GDMATE accessible to
coding novices. For those of you who have little experience with Python,
here are some specific features and pitfalls of the language:

* Python uses specific indentation. A script might fail if a code block is not indented correctly. We use four spaces and no tabs. Mixing spaces and tabs can cause trouble.
* Indices should be given inside square brackets and function or method call arguments inside parentheses (different from Matlab).
* The first index of an array or list is 0 (e.g. x[0]), not 1.
* Put dots after numbers to make them floats instead of integers.

## Contributing to GDMATE
GDMATE is a community-driven, open-source Python package by and for the geodynamics community. If you have code you would like to contribute, please review our [contribution guidelines](https://gdmate.readthedocs.io/en/latest/CONTRIBUTING.html) and open a [pull request](https://docs.github.com/en/pull-requests).
