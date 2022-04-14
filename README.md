# countland

Tools for analyzing biological count data, especially from single cell RNA-seq. 

Please see our manuscript for more details:

> Church et al. 2022. Counting RNA transcripts in cells. XXXX


`countland` is implemented in both R and python. The code for each is included in this repository.

## python

### Installation for python

To prepare a conda environment and install countland (before first use):

    conda create -n countland -c conda-forge
    conda activate countland
    pip install git+https://github.com/shchurch/countland.git#egg=countland&subdirectory=countland-py

To activate the conda environment (before each use):

    conda activate phylopytho

### Development in python

To install from local source:

    conda create -n countland -c conda-forge python==3.10
    conda activate countland
    cd countland-py
    pip install .[dev]
    pip install jupyter

For testing, this will reload countland functions (e.g. after an edit)

    import importlib
    importlib.reload(countland)

#### Running tests in python

In the `countland-py/`, run:

    conda activate countland
    python -m pytest


## R

### Installation for R

From an R prompt, run the following: 

    library(devtools)
    install_github("shchurch/countland", subdir="countland-R")


### Development in R

This package is built with the excellent [devtools](https://github.com/hadley/devtools). Extensive explanations on using devtools
to create R packages is available in the book
[R Packages](http://r-pkgs.had.co.nz/).

Development typically involves `cd`ing to the package directory `countland-R/`, launching R, and running some combination of the following:

	options(error=traceback) # Get line numbers for errors
    library(devtools)
    load_all()
    test()
    document()
    check() # A wrapper for R CMD check, see http://r-pkgs.had.co.nz/check.html#check
    build() # Create package bundle, including executed vignettes

To regenerate the pdf manual, run the following shell command in the package directory:

    R CMD Rd2pdf . --force --output=countland-manual.pdf