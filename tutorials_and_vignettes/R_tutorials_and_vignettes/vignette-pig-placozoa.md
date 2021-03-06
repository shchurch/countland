countland: pig and placozoa
================

### Samuel H. Church

This document applies `countland` to the analysis of real-world
scRNA-seq data.

The following packages are required to run the analysis:

``` r
library(countland)
library(Seurat)
library(aricode)
library(clevr)

library(knitr)
library(dplyr)
library(gridExtra)

library(ggplot2)
library(viridis)
color_palette <- c("#8c564b", "#9467bd", "#2ca02c", "#e377c2", "#d62728", "#17becf", "#bcbd22", "#ff7f0e", "#7f7f7f", "#1f77b4")
theme_set(theme_classic())
```

## Get the data

Here we apply countland to a dataset from the placozoan (*Trichoplax
adherens*) provided by [Sebé-Pedros *et al*
(2018)](https://doi.org/10.1038/s41559-018-0575-6), and a dataset from
the pig (*Sus scrofa*) provided by [van Zyl *et al*
(2020](https://www.pnas.org/cgi/doi/10.1073/pnas.2001250117).

``` r
Tri.data <- Seurat::Read10X(data.dir = "../data/Trichoplax_SebePedros2018/")
Sus.data <- Seurat::Read10X(data.dir = "../data/Pig9_vanZyl2020_mtx/")
```

## Run countland

We first apply spectral embedding to the two datasets.

``` r
Sus <- countland(Sus.data,remove_empty = TRUE)
Sus <- Subsample(Sus,gene_counts = ncol(Sus@counts)*50,cell_counts = "min")
Sus <- Dot(Sus,subsample=TRUE)
Sus <- Embed(Sus,n_components=20)

Tri <- countland(Tri.data,remove_empty = TRUE)
Tri <- Subsample(Tri,gene_counts = ncol(Tri@counts)*50,cell_counts = "min")
Tri <- Dot(Tri,subsample=TRUE)
Tri <- Embed(Tri,n_components=20)
```

Because we have no *a priori* hypothesis as to the number of clusters
that may be present, we check to see if the eigengap heuristic for
spectral clustering is informative.

``` r
PlotEigengap(Sus)
```

![](vignette-pig-placozoa_files/figure-gfm/elbow-Sus-1.png)<!-- -->

For *S. scrofa*, it appears as though the largest visible gap at the
highest value is between 7 and 8, so we’ll use 7 clusters.

``` r
PlotEigengap(Tri)
```

![](vignette-pig-placozoa_files/figure-gfm/elbow-Tri-1.png)<!-- -->

For *T. adherens*, it appears that the largest visible gap is between 11
and 12, so we’ll use 11 clusters.

``` r
Sus <- Cluster(Sus,n_clusters=7,n_components=7)
#> [1] "Performing spectral clustering on dot products..."
#> [1] "    done."
Tri <- Cluster(Tri,n_clusters=11,n_components=11)
#> [1] "Performing spectral clustering on dot products..."
#> [1] "    done."
```

We can visualize the results using spectral embedding.

``` r
embed <- Tri@embedding[,2:3]
embed <- setNames(data.frame(embed),paste("component_",seq_len(2),sep=""))

ggplot2::ggplot(embed,aes(x = component_1,y = component_2, color=as.character(Tri@cluster_labels))) +
geom_point(size=1) +
guides(color=guide_legend(title="cluster")) +
scale_color_manual(values=c(color_palette,"black"))
```

<img src="vignette-pig-placozoa_files/figure-gfm/plot-Tri-1.png" style="display: block; margin: auto;" />

In the case of *S. scrofa* we see that the first two variable components
from sepctral embedding are not very informative, so it’s helpful to
view additional components.

``` r
Scl <- paste0("countland_cluster:",as.character(Sus@cluster_labels))
Scolor <- viridis::turbo(n=8)
Sdf <- setNames(data.frame(Sus@embedding[,2:5]),paste0("component_",seq_len(4)))
Sdf$clusters <- Scl
g1 <- ggplot(Sdf,aes(x = component_1, y = component_2, color = clusters)) + 
  geom_point(size=1) + 
  scale_color_manual(values=Scolor) +
  theme("legend.position" = "none") + 
  theme(axis.ticks = element_blank(),axis.text = element_blank())
g2 <- ggplot(Sdf,aes(x = component_3, y = component_4, color = clusters)) + 
  geom_point(size=1) + 
  scale_color_manual(values=Scolor) +
  theme("legend.position" = "none") + 
  theme(axis.ticks = element_blank(),axis.text = element_blank())
grid.arrange(g1,g2,ncol=2)
```

<img src="vignette-pig-placozoa_files/figure-gfm/plot-Sus-1.png" style="display: block; margin: auto;" />
