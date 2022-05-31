# countland

Tools for analyzing biological count data, especially from single cell RNA-seq. 

Please see our manuscript for more details:

> Church et al. 2022. Normalizing need not be the norm: count-based math for analyzing single-cell data. _in prep_

`countland` is implemented in both `R` and `python`. The code for each is included in this repository.

## `python`

### Installation for `python`

To prepare a `conda` environment and install `countland` (before first use):

    conda create -n countland -c conda-forge
    conda activate countland
    pip install git+https://github.com/shchurch/countland.git#subdirectory=countland-py

To activate the conda environment (before each use):

    conda activate countland

### Running the tutorial in `python`

The easiest way to run the tutorial is as a Google Colab notebook. Just open the following link and follow the instructions:

https://colab.research.google.com/github/shchurch/countland/blob/main/countland-py/vignettes/vignette_tutorial.ipynb

Alternatively, [the `python` tutorial](./countland-py/vignettes/vignette_tutorial.ipynb) can be run locally in a [jupyter notebook](https://jupyter.org/).

## `R`

### Installation for `R`

From an `R` prompt, run the following: 

    library(devtools)
    install_github("shchurch/countland", subdir="countland-R")

### Running the tutorial in `R`

[The `R` tutorial](./countland-R/vignettes/vignette-tutorial.Rmd) can be run locally as an Rmarkdown file, e.g. knit in [RStudio](https://www.rstudio.com/).

## Development

See [development.md](./development.md) for details.
