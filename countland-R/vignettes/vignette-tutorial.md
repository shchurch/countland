countland: R tutorial
================

### Samuel H. Church

This tutorial demonstrates the major functions of `countland` by
applying them to a Gold standard single-cell RNA sequencing dataset from
[Freytag *et al* (2018)](10.12688/f1000research.15809.2).

`countland` is a **barebones** set of functions for applying a
restricted linear algebra to the analysis of count-based data. As such,
there are many opportunities for further optimization that may prove
useful in the anlaysis of your own data. We provide the source code
freely available at <https://github.com/shchurch/countland> and
encourage users and developers to fork the code for their own purposes.

The following packages are required to complete the tutorial

``` r
library(countland)
library(Seurat)
#> Warning: package 'Seurat' was built under R version 4.1.2
#> Attaching SeuratObject
#> Attaching sp
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.1.2
theme_set(theme_classic())
```

## Get the data

We have used the Gold standard dataset provided by [Freytag *et al*
(2018)](10.12688/f1000research.15809.2). This data consists of \~1000
cells that have ground truth labels corresponding to three human lung
cell lines.

`countland` accepts an sparse data matrix (as does `Seurat`).

``` r
gold_path <- system.file("exdata", "Gold_Freytag2018/", package = "countland", mustWork = TRUE)
gold.data <- Seurat::Read10X(data.dir = gold_path)
m <- gold.data
```

## Initialize `countland` object

Initialize `countland` by calling the core function on sparse data
matrix.

``` r
C <- countland(m)
#> [1] "countland object"
#> [1] "after removing empty cells and genes,"
#> [1] "the count matrix has 29212 genes (rows)"
#> [1] "    and 925 cells (columns)"
#> [1] "the fraction of entries that are nonzero is 0.2905"
```

The count matrix is stored in `C@counts`

``` r
C@counts[1:10,1:10]
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>    [[ suppressing 10 column names 'H2228_AAACCTGCAGACACTT-1', 'H1975_AAAGATGCACATTTCT-1', 'H1975_AAAGATGTCCTTTACA-1' ... ]]
#>                                    
#> ENSG00000243485 . . . . . . . . . .
#> ENSG00000238009 . . . . . . . . . .
#> ENSG00000239945 . . . . . . . . . .
#> ENSG00000233750 . . . . . . . . . .
#> ENSG00000268903 . . . . . . . . . .
#> ENSG00000241860 . . . . . . . 1 . .
#> ENSG00000279457 . . . . . . . . . .
#> ENSG00000228463 . . . . . 1 . . . .
#> ENSG00000236601 . . . . . . . . . .
#> ENSG00000237094 . . . . . 2 . . . .
```

Note that most counts are zero for scRNA-seq data.

## Summarize counts across cells and genes

We can explore our data by checking the total number of counts and other
expression measures across cells.

Expression measures include: \* total counts \* maximum number of counts
observed across genes \* number of genes with observed counts \* number
of genes with counts above 1, or 10 \* number of unique count values per
cell \* count index = number of *n* genes ≥*n* counts

It can be helpful to see how many counts are derived from certain genes,
such as mitochondrial genes (here we have chosen a gene name with
particularly high counts, regular expressions are permitted for name
matching).

``` r
# sum counts from genes matching this regex expression
gene_string <- "162396$"

# calculate expression scores across cells
C <- ScoreCells(C,gene_string)
head(C@cell_scores)
#>                      names n_features max_count_value total_counts
#> 1 H1975_AAAGATGCACATTTCT-1       8574            1010       100763
#> 2 H1975_AAAGATGTCCTTTACA-1       9190            1803       121201
#> 3 H1975_AAATGCCCACTTCGAA-1       8433            1097        95220
#> 4 H1975_AAATGCCTCATATCGG-1       8311            1678        99321
#> 5 H1975_AACACGTGTCAGAAGC-1       8642            1142       106598
#> 6 H1975_AACCATGAGCGTAGTG-1       7907            1549        90917
#>   n_features_above1 n_features_above10 unique_count_values count_index
#> 1              5746               1432                 227         127
#> 2              6316               1567                 245         137
#> 3              5643               1390                 219         127
#> 4              5556               1247                 221         131
#> 5              5833               1431                 231         133
#> 6              5058               1127                 214         128
#>   feature_match_counts
#> 1                    0
#> 2                    1
#> 3                    0
#> 4                    0
#> 5                    0
#> 6                    0
```

We can also calculate the same meaures, but across cells.

``` r
# calculate expression scores across genes
C <- ScoreGenes(C)
head(C@gene_scores)
#>             names n_cells max_count_value total_counts n_cells_above1
#> 1 ENSG00000000003     825              26         3289            662
#> 2 ENSG00000000419     922              58        10379            914
#> 3 ENSG00000000457     226               4          278             40
#> 4 ENSG00000000460     446               6          805            202
#> 5 ENSG00000000938       2               1            2              0
#> 6 ENSG00000000971     321               7          565            144
#>   n_cells_above10 unique_count_values count_index
#> 1              38                  20          13
#> 2             384                  43          30
#> 3               0                   4           3
#> 4               0                   6           5
#> 5               0                   1           1
#> 6               0                   7           6
```

## Cluster cells by similarity

The dot (or inner) product is a measure of similarity between vectors.
In this case, it tells us how similar two cells are based on the
distribution of transcript counts, and scaled by the total counts per
cell. A dot product of 0 indicates orthogonal cell vectors (no shared
counts), larger values indicate aligned cell vectors.

``` r
C <- Dot(C)
#> [1] "Calculating dot products between rows..."
#> [1] "    done."
```

Cell populations can be compared and distinguished by embedding and
clustering the matrix of pairwise dot products (contained in `C.dots`).
This matrix is an unbounded affinity matrix. It is symmetric, and
contains only integer values above 0. Spectral embedding clustering is
appropriate for this type of matrix.

First, we embed to investigate the optimal number of clusters for our
data.

``` r
C <- Embed(C)
#> [1] "Performing spectral embedding on dot products..."
#> [1] "    done."
```

The eigengap heuristic is can help decide the optimal number of
clusters, but it is only a guideline. According to this heuristic, the
optimal number of clusters is *k* where the difference in eigenvalues $
\| e\_{k+1} - e\_{k} \| $ is largest.

``` r
PlotEigengap(C)
```

![](vignette-tutorial_files/figure-gfm/eigengap-1.png)<!-- -->

For many datasets, you may want to consider other factors, e.g. choosing
a minimum number of clusters, whether or not the eigengap reflects this.

Here we have chosen 3 as the optimal number of clusters.

``` r
C <- Cluster(C,n_clusters=3)
#> [1] "Performing spectral clustering on dot products..."
#> [1] "    done."
```

We can now visualize clusters using spectral embedding.

``` r
PlotEmbedding(C)
```

![](vignette-tutorial_files/figure-gfm/plot-cluster-1.png)<!-- -->

## Subsampling data

Cells are not sequenced to standard sequencing depth. This is sometimes
a problem for downstream comparisons, but not always.

You can create an alternative count matrix where cells have an equal
number of counts using the function `Subsample()`. To subsample to the
minimum number of total counts across cells, use `cell_counts='min'`,
otherwise the number of counts must be larger than the minimum value.

This matrix is stored in `C@subsample`

``` r
C <- Subsample(C,cell_counts='min')
#> [1] "subsampling all cells to a standard sequencing depth of 40152"
```

Similarly, there is often substantial heterogeneity in the magnitude of
expression across genes. This may result in highly expressed genes
having an outsized impact on results.

You can create an alternate count matrix where gene expression is
bounded at a maximum total counts across cells, use
`Subsample(gene_counts=[maximum]`. An example maximum value might be
equal to 10x the number of cells (columns in the count matrix).

``` r
C <- Subsample(C,gene_counts=10*ncol(C@counts)) # doing so will overwrite our previous subsampled matrix 
#> [1] "subsampling 1592 genes to a max total counts of 9250"

# to subsample both genes and cells, use gene_counts and cell_counts in the same function. Genes will be subsampled first.

#C <- Subsample(C,gene_counts=10*ncol(C@counts),cell_counts='min')
```

## Subsetting data

You can filter the count matrix to only certain cells and genes using
`SubsetCells()` and `SubsetGenes()`.

``` r
filter_cell_names <- C@cell_scores[C@cell_scores$n_features < 8500,]$names
filter_cell_index <- which(C@names_cells %in% filter_cell_names) # cells with fewer than 8,500 unique features
C <- SubsetCells(C,filter_cell_index,remove_empty=FALSE)
#> [1] "new number of genes: 29212"
#> [1] "new number of cells: 463"
```

``` r
filter_gene_names <- C@gene_scores[C@gene_scores$n_cells > 100,]$names
filter_gene_index <- which(C@names_genes %in% filter_gene_names) # cells with greater than 100 observations
C <- SubsetGenes(C,filter_gene_index,remove_empty=FALSE)
#> [1] "new number of genes: 13368"
#> [1] "new number of cells: 463"
```

With `countland`, such data filtering may not be necessary or helpful.
The original count matrix can be restored at any time.

``` r
C <- RestoreCounts(C)
#> [1] "new number of genes: 29212"
#> [1] "new number of cells: 925"
```

## Identify marker genes

What makes a gene an ideal marker for a cluster may depend on downstream
applications. For example, the ideal marker gene might be defined as the
gene detected in all cells in a given cluster and none of the rest.
Under this definition, the top marker gene for each cluster can be
identified by counting and comparing the number of non-zero cells.

``` r
C <- RankMarkerGenes(C,method='prop-zero',subsample=F)
C@marker_genes[(C@marker_genes$cluster_label == 1),]
#>                           names gene_index diff.zeros rank cluster_label
#> ENSG00000253706 ENSG00000253706      13393  0.9155633    1             1
#> ENSG00000258484 ENSG00000258484      21428  0.8744115    2             1
#> ENSG00000233429 ENSG00000233429      10683  0.8721458    3             1
#> ENSG00000167641 ENSG00000167641      27283  0.8431677    4             1
#> ENSG00000225548 ENSG00000225548       4990  0.8379180    5             1
#> ENSG00000169213 ENSG00000169213        933  0.8241418    6             1
#> ENSG00000260027 ENSG00000260027      24360  0.8069948    7             1
#> ENSG00000249395 ENSG00000249395      13397  0.7898422    8             1
#> ENSG00000134709 ENSG00000134709       1030  0.7687606    9             1
#> ENSG00000261780 ENSG00000261780      25457  0.7581397   10             1
```

``` r
cluster_marker <- C@marker_genes[C@marker_genes$cluster_label == 1,]
cluster_top <- cluster_marker[cluster_marker$rank == 1,]$names
gene_index = which(C@names_genes == cluster_top)
PlotMarker(C,gene_index)
```

![](vignette-tutorial_files/figure-gfm/plot-marker-zero-1.png)<!-- -->

Alternatively, the top marker genes for each cluster can be identified
by ranking genes according to differential gene expression, calculated
using the Wilcoxon rank-sum statistic.

When calculating differential gene expression, it typically makes sense
to subsample cells to a standard sequencing depth.

``` r
C <- Subsample(C,cell_counts='min')
#> [1] "subsampling all cells to a standard sequencing depth of 40152"
C <- RankMarkerGenes(C,method='rank-sums',subsample=T)
C@marker_genes[(C@marker_genes$cluster_label == 1),]
#>                 names gene_index statistic        pvalue    adj.pvalue
#> 27283 ENSG00000167641      27283  177633.0 4.733655e-168 1.382795e-163
#> 1550  ENSG00000232527       1550  172611.0 2.449540e-138 3.120273e-134
#> 18498 ENSG00000166908      18498  180148.5 3.204443e-138 3.120273e-134
#> 18560 ENSG00000183735      18560  178794.0 1.072235e-133 7.830535e-130
#> 11545 ENSG00000128591      11545  170106.0 9.733681e-133 5.686806e-129
#> 14061 ENSG00000147889      14061  180848.0 5.965072e-132 2.904195e-128
#> 10909 ENSG00000132434      10909  177715.0 3.371481e-131 1.406967e-127
#> 18504 ENSG00000135506      18504  180579.5 2.129375e-130 7.775414e-127
#> 10902 ENSG00000132432      10902  180957.0 7.726720e-130 2.436402e-126
#> 10910 ENSG00000154978      10910  180405.5 8.340414e-130 2.436402e-126
#>       significant rank cluster_label
#> 27283        TRUE    1             1
#> 1550         TRUE    2             1
#> 18498        TRUE    3             1
#> 18560        TRUE    4             1
#> 11545        TRUE    5             1
#> 14061        TRUE    6             1
#> 10909        TRUE    7             1
#> 18504        TRUE    8             1
#> 10902        TRUE    9             1
#> 10910        TRUE   10             1
```

``` r
cluster_marker <- C@marker_genes[C@marker_genes$cluster_label == 1,]
cluster_top <- cluster_marker[cluster_marker$rank == 1,]$names
gene_index = which(C@names_genes == cluster_top)
PlotMarker(C,gene_index)
```

![](vignette-tutorial_files/figure-gfm/plot-marker-ranks-1.png)<!-- -->
