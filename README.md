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



## R

### Installation for R

From an R prompt, run the following: 

    library(devtools)
    install_github("shchurch/countland", subdir="countland-R")

## Development

See [development.md](./development.md) for details on development.