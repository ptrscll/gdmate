# Building a Python Package #

As with all things Python, there are many ways to create a Python package. There are guides to many of the individual tools and steps on the Internet, but it can be difficult to find a clear path from zero to package, so we will attempt to document how GDMATE was created here.

This guide assumes basic familiarity with Python and installing packages.

## Basic directory structure ##

The root directory contains the usual basic repository files (`README.md`, `license.txt`, `.gitignore`), as well as the `setup.py` file for the package and some additional setup files described below. The actual source code for the installable package resides in the `gdmate` directory, which has the same name as the repository so that the package will also be named `gdmate`. Additional directories contain material that isn't part of the installable source code, including documentation (`docs`), Jupyter Notebooks (`notebooks`), and tests (`tests`), which are each discussed below.

Within the `gdmate` directory, the source code is contained within _modules_, which are individual `.py` files that each contain callable _functions_. These modules are organized into _packages_, which are directories containing multiple modules and an `__init__.py` file, which indicates that the modules should be treated together as a package. The principal package is `gdmate`, and each of the subdirectories (e.g. `analysis_modules`, `io`, etc.) are considered a _subpackage_ of gdmate. Each package/subpackage needs to be installed during package setup, and the package/subpackage structure is an important consideration for importing and namespaces wihen using the package.

More about Python modules available [here](https://docs.python.org/3/tutorial/modules.html), and more about packaging [here](https://packaging.python.org/en/latest/tutorials/packaging-projects/).

## Setup file for pip installation ##

At its core, Python package installation is usually handled using pip (even in conda environments). When the command `pip install` is executed, pip will look for a setup file, with one of the standard options being a `setup.py` file that is configured using the setuptools package. A sample file from the Python Packaging Authority is available [here](https://github.com/pypa/sampleproject/blob/main/setup.py).

Fundamentally, this file contains key metadata (authors, version, etc.), information about the location of source code, and information about dependencies. In GDMATE, the following lines are key for setting up the package:

```
    packages=["gdmate","gdmate.visualization","gdmate.analysis_modules",
        "gdmate.io","gdmate.education","gdmate.material_models"],
    python_requires=">=3.7, <4",
    install_requires=["numpy","scipy","matplotlib","pyvista"],
```

The `packages` variable specifies where to find the source code and includes both the main package in the main directory and subpackages within subdirectories.

The `python_requires` variable specifies the supported version of Python. In this case, we are trying to maintain compatibility with all version of Python 3 from 3.7 on, which is reflected in our testing workflow discussed below.

The `install_requires` variable lists the dependencies for the package. Note that in many cases we can use the dependencies specified here in place of a `requirements.txt` file.

With this file in place, running `pip install .` in the root directory will install GDMATE in the current Python environment, because pip will search that directory and use the information in `setup.py` to install the package.

## Imports and Namespaces ##
If all `__init__.py` files are blank, each subpackage within the Python package can be imported directly (e.g., `import gdmate.analysis_modules`). However, if only the root package is imported (i.e., `import gdmate`), the modules within subpackages will not be accessible. Adding import statements to the base `__init__.py` file defines the namespaces for subpackages, modules, and or functions in relation to the base package. For example, the `__init__.py` file for GDMATE contains the line:

```
from gdmate.visualization import pyvista_vis
```

As a result, when a user's script runs `import gdmate`, the `pyvista_vis` module can be accessed simply as `gdmate.pyvista_vis`. Note that this particular formulation removes the subpackage `visualization` from the namespace; `gdmate.visualization.pyvista_vis` will fail unless the import statement is `import gdmate.visualization`. Python namespaces are notoriously confusing.


## Testing ##
Setting up automated tests is essential for debugging non-functional code and ensuring compatibility with multiple versions of Python. We have a relatively simple testing workflow implemented using `pytest` and automated using GitHub Actions.

### Pytest ###
Designing tests for use with Pytest is fairly straightforward. Pytest will search a repository for directories, modules, and functions with the word "test," making it simple to run all tests just with the command `pytest`. To design a test for use with Pytest, you simply have to make functions with `assert` statements that Pytest can attempt to evaluate as true. A simple example is shown [here](https://docs.pytest.org/en/7.1.x/). In GDMATE, the tests are housed within a `tests` directory separate from the source code, and the tests can be run locally by installing Pytest in an environment with GDMATE and executing `pytest`.

We extend our testing to include whether all code in Jupyter Notebooks will also execute successfully, which is very easy to implement with the nbmake plugin for Pytest. Once nbmake is installed in the environment, all that is required to include Jupyter Notebooks is to run `pytest --nbmake`. The Jupyter Notebooks for GDMATE are all housed in the `notebooks` directory outside of the source code.

Although tests can be run locally, we automatically run tests using multiple versions of Python a GitHub Actions workflow. This workflow is contained within `.github/workflows/python-test.yml` and contains instructions to set up a virtual environment for versions 3.7, 3.8, 3.9, and 3.10 of Python, install GDMATE along with its dependencies, and Pytest, and then run Pytest. These tests are run every time there is a commit or pull request on the `main` branch of the GDMATE repository to ensure that new additions do not break the package.

## Sphinx Documentation ##
Documentation is generated using Sphinx and hosted online using Read the Docs. The necessary files are housed in a `docs` directory separate from the source code. As described [here](https://www.sphinx-doc.org/en/master/tutorial/getting-started.html), the essential files needed for Sphinx can be generated using `sphinx-quickstart` once Sphinx is installed in an environment. The most important of these is `conf.py`, which contains the information needed to build the documentation, including setting the package directory, choosing the theme, and specifying any extensions needed. We make use of the `autosummary` and `autodoc` extensions included within Sphinx to autodocument the package API, as well as `nbsphinx` for rendering Jupyter Notebooks and `myst-parser` for rendering Markdown files.
