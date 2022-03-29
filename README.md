# countland
tools in python and R for analyzing biological count data, especially from single cell RNAseq 


## Installation for R

From an R prompt, run the following: 

    install_github("shchurch/countland", subdir="countland-R", auth_token="TOKEN")

Replace `TOKEN` with a [github token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token). (This can be removed once the repo is public)


## Installation for python

To prepare a conda environment and install countland (before first use):

    conda create -n countland -c conda-forge
    conda activate countland
    pip install git+https://github.com/shchurch/countland.git#egg=countland&subdirectory=countland-py

To activate the conda environment (before each use):

    conda activate phylopytho

## Development in R

## Development in python

To install from local source:
    conda create -n countland -c conda-forge
    conda activate countland
    cd countland-py
    pip install .[dev]

For testing, this will reload countland functions (e.g. after an edit)

    import importlib
    importlib.reload(countland)

### Running tests

In the `countland-py/`, run:

    conda activate countland
    python -m pytest