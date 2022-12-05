# Contributing Code to GDMATE #
The underlying goal of GDMATE is to enable highly transparent, community-drive code development within a structure that keeps the code organized, accessible, and useable. With that in mind, we welcome all contributions and ask that code abide by the following guidelines. Please contact one of the maintainers if you have questions about how to organize, document, and test your code for inclusion in the package. These guidelines will be used in assessing all pull requests for new code.

## Summary Checklist for Code Contributions ##
* Code written as functions in modules within subpackages.
    - `def` statements in `.py` files
* Documentation in docstring, comments, and Sphinx documentation.
    - docstring enclosed in 2 sets of `"""` after each `def` statement
    - comments preceded by `#` before lines of code
    - modules added to `docs/api.rst`
*  For new modules, `import` statement added to `gdmate/__init__.py`.
* Test function for each function added to corresponding test module in `tests`.
* Example code usage added and annotated in Jupyter Notebook under `notebooks`.


## Functional-Style Python Modules ##
Code contributions should be organized into Python _functions_ (which start with a `def` statement) organized into _modules_ (`.py` files). Each module should be placed within a _subpackage_, one of the first-order directories within the `gdmate` source code directory. Contributions may include adding functions to existing modules or creating new modules.

## Documentation ##
All new functions should contain a _docstring_, which appears in the code immediately following the `def` statement and is encased in two sets of three double quotation marks (`"""`). The docstring should contain a one-line description of the function, followed by a separate paragraph with any additional details and then a list of all parameters and returns. Below is an example docstring, and additional information can be found [here](https://peps.python.org/pep-0257/).
```
def function(A,B):
    """
    This is what the function does.

    This is more information about what the function does.

    Parameters:
        A : type
            This is what A is.
        B : type
            This is what B is.

    Returns:
        X : type
            This is what X is.
    """
```
To add the information from the docstring to the Sphinx documentation for GDMATE, the module in which the function appears must appear in the file `docs/api.rst`. If you are adding functions to an existing module, this should be handled automatically. However, new modules will need to be added as a line to this file:

```
.. autosummary::
    :toctree: generated 
    :template: module.rst
    :recursive:

    gdmate.module1
    gdmate.module2
```

Additionally, lines of code where the functionality of the code is not immediately apparent should be preceded by _comments_ describing what the code does. Comments are preceded by the `#` symbol:
```
# This converts the file to a mesh object.
mesh = pv.read(file)
```
## The `__init__.py` file ##
The current structure of GDMATE is such that users should be able to call each module after executing `import gdmate`. For this to work, an import statement for each new module needs to appear in the `__init__.py` file in the main `gdmate` source code directory, using the following format:
```
from gdmate.subpackage import module1
from gdmate.subpackage import module2
```
This will allow functions within the module to be called in the following way within a script:
```
import gdmate as gd

gd.module1.function()
```

## Tests ##
Every function within a module needs a corresponding test to ensure that it is operating as intended when changes are made to the code and across Python versions. Each module should have a `.py` file in the `tests` directory with the name after the format `test_module1.py`. Within the file, a test function needs to be defined for each function in the module. The test function should call the function and check that the output is as intended using `assert` statements:
```
from gdmate.subpackage import module1

def test_function():
    output = module1.function(input)
    assert output == correct_answer
```
## Jupyter Notebooks ##
It is also important that the functionality of each module be documented in a Jupyter Notebook so users can quickly see examples of how to use GDMATE. Notebooks do not necessarily need to be organized exactly according to the subpackage/module/function hierarchy, but all new code should have corresponding annotated content in a new or existing Jupyter Notebook (in the `notebooks` directory) that demonstrates the key features.
