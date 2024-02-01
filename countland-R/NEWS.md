# countland 0.1.2

* `ScikitManifoldSpectralEmbedding()`, the spectral embedding utility function, now works with behavior of `igraph` v. 2.0.0, by removing self-loops in the adjacency matrix (setting `diag = F`). In addition the function was changed from the deprecated `igraph::graph.adjacency()` to `igraph::graph_from_adjacency_matrix()` (#18).

# countland 0.1.1

* Added a `NEWS.md` file to track changes to the package.

