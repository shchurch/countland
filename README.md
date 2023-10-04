# countland

Tools for analyzing biological count data, especially from single cell RNA-seq. 

Please see our [manuscript](https://www.biorxiv.org/content/10.1101/2022.06.01.494334v1) for more details:

> Church et al. 2022. Normalizing need not be the norm: count-based math for analyzing single-cell data. _bioRxiv_ https://doi.org/10.1101/2022.06.01.494334

`countland` is implemented in both `R` and `python`. The code for each is included in this repository.

## `python`

### Installation for `python`

`countland` is available with `pip`: https://pypi.org/project/countland/

To prepare a `conda` environment and install `countland` (before first use):

    conda create -n countland -c conda-forge
    conda activate countland
    pip install countland

To activate the conda environment (before each use):

    conda activate countland

The develompent version from in this repository can be installed using 

    pip install git+https://github.com/shchurch/countland.git#subdirectory=countland-py

### Running the tutorial in `python`

The easiest way to run the tutorial is as a Google Colab notebook. Just open the following link and follow the instructions:

https://colab.research.google.com/github/shchurch/countland/blob/main/tutorials_and_vignettes/python_tutorials_and_vignettes/vignette_tutorial.ipynb

Alternatively, [the `python` tutorial](./tutorials_and_vignettes/python_tutorials_and_vignettes//vignette_tutorial.ipynb) can be run locally in a [jupyter notebook](https://jupyter.org/).

## `R`

### Installation for `R`

`countland` is available from CRAN: https://CRAN.R-project.org/package=countland

From an `R` prompt, run the following: 

    install.packages("countland")

Teh development version from this repository can be installed using 

    library(devtools)
    install_github("shchurch/countland", subdir="countland-R")

### Running the tutorial in `R`

[The `R` tutorial](./tutorials_and_vignettes/R_tutorials_and_vignettes//vignette-tutorial.Rmd) can be run locally as an Rmarkdown file, e.g. knit in [RStudio](https://www.rstudio.com/).

## Development

See [development.md](./development.md) for details.
