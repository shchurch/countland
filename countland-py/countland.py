import sys
import logging

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

from sklearn.cluster import SpectralClustering
from sklearn.manifold import SpectralEmbedding
from sklearn.metrics import pairwise
from sklearn import preprocessing

from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu
from statsmodels.stats import multitest

from integer_matrix_approximation import *

class countland:
    """
    A class used to represent a countland RNA count object
    ...

    Attributes
    ----------
    counts : numpy.ndarray
        count matrix with rows as genes, columns as cells
    names_cells : numpy.ndarray
        cell names
    names_genes : numpy.ndarray
        gene names
    raw_counts : numpy.ndarray
        count matrix as originally loaded
    raw_names_cells : numpy.ndarray
        cell names as originally loaded
    raw_names_genes : numpy.ndarray
        gene names as originally loaded

    Methods
    -------
    RestoreCounts()
        Restores count matrix to original numpy array
    SubsetGenes(gene_indices)
        Keeps only genes listed in numpy array
    SubsetGenes(gene_indices)
        Keeps only cells listed in numpy array
    Subsample(n_counts)
        Subsamples cells to standard number of total counts
    ScoreGenes(subsample=True)
        Calculates several expression scores based on counts
    Dot()
        Calculates pairwise dot products between all cells based on counts
    Cluster(n_clusters)
        Performs spectral clustering using dot products between cells
    PlotEmbedding()
        Plots cell clusters using spectral embedding
    RankMarkerGenes(method='prop-zero')
        Ranks marker genes for each cluster based on expression metrics
    PlotMarker(gene_index)
        Plots marker gene expression across cells using spectral embedding
    PlotIMAElbow(max_features,u_bounds,subsample=True)
        Approximates count matrix using a range of feature numbers and plots results
    IMA(features,u_bounds,l_bounds=[0,0],maxiter=1000000.0,stop_crit=0.0001,subsample=True)
        Uses integer matrix approximation to identify matrices with a set number of features that can be multiplied to approximate the count matrix
    PlotIMA(x=0,y=1,subsample=True)
        Plots cells using two features from integer matrix approximation
    SharedCounts(n_clusters,n_cells=100,subsample=True)
        Clusters and collapses genes based on number of shared counts
    PlotIMA(x=0,y=1,subsample=True)
        Plots cells using two features from the collapsed gene matrix
    """
    
    def validate(self):
        assert self.counts.shape[0] == len(self.names_cells), "number of rows differs from number of cell names"
        assert self.counts.shape[1] == len(self.names_genes), "number of columns differs from number of gene names"
        return True
    
    def __init__(self,a,verbose=True):
        if verbose == True:
            logger = logging.getLogger()
            logger.handlers = []
            logger.setLevel(logging.INFO)
            
            handler = logging.StreamHandler(sys.stdout)
            
            logger.addHandler(handler)
        else:
            logger = logging.getLogger()
            logger.handlers = []
            logger.setLevel(logging.WARNING)
            
            handler = logging.StreamHandler(sys.stdout)
            
            logger.addHandler(handler)            
            
        logging.info("Initializing countland object")
        # parse the AnnData object
        self.counts = a.X.toarray().astype(int)
        assert np.linalg.norm(a.X.toarray() - self.counts) < 0.00001, "non-integer counts present"
        
        self.names_cells = np.array(a.obs_names)
        self.names_genes = np.array(a.var_names)
        
        self.raw_counts = np.copy(self.counts)
        self.raw_names_cells = np.copy(self.names_cells)
        self.raw_names_genes = np.copy(self.names_genes)

        self.validate()        
 
    def __str__(self):
        return f'''
        countland object
        Count matrix has {self.counts.shape[0]} cells (rows) 
         and {self.counts.shape[1]} genes (columns)
        The fraction of entries that are nonzero is {round(self._FractionNonzero(),4)}
        '''
    
    def _FractionNonzero(self):
        """
        Calculates fraction of non-zero counts
        """ 
        n = self.counts.shape[0] * self.counts.shape[1]
        return np.count_nonzero( self.counts ) / n  
              
    def RestoreCounts(self):
        """
        Restores count matrix to original state
        ...
        
        Updates
        ----------
        counts : numpy.ndarray
        """ 
        self.counts = np.copy(self.raw_counts)
        self.names_cells = np.copy(self.raw_names_cells)
        self.names_genes = np.copy(self.raw_names_genes)
        
        logging.info("Number of genes: %s",self.counts.shape[1])
        logging.info("Number of cells: %s",self.counts.shape[0])
        
    def SubsetGenes(self,gene_indices):
        """
        Subsets genes using a numpy array of gene indices  
        ...
        
        Parameters
        ----------
        gene_indices : numpy.ndarray
            index values of genes to keep 
    
        
        Updates
        ----------
        counts : numpy.ndarray
        """ 
        
        assert isinstance(gene_indices,np.ndarray), "expecting numpy array"
        self.counts = self.counts[:,gene_indices]
        self.names_genes = np.array(self.names_genes)[gene_indices]
        
        logging.info("New number of genes: %s",self.counts.shape[1])

    def SubsetCells(self,cell_indices):
        """
        Subsets cells using a numpy array of cell indices  
        ...
        
        Parameters
        ----------
        gene_indices : numpy.ndarray
            index values of cells to keep 
    
        
        Updates
        ----------
        counts : numpy.ndarray
        """ 
        
        assert isinstance(cell_indices,np.ndarray), "expecting numpy array"
        self.counts = self.counts[cell_indices,:]
        self.names_cells = np.array(self.names_cells)[cell_indices]
        
        logging.info("New number of cells: %s",self.counts.shape[0])

    def _CountIndex(self,c):
        """
        Internal function for calculating count index
        """ 
        
        unique,counts = np.unique(c,return_counts=True)
        return(np.max(unique[counts >= unique]))

    def ScoreCells(self,gene_string=None):
        """
        Calculates several scores for counts across cells
        ...
        
        Parameters
        ----------
        gene_string : str
            Regular expression pattern matching gene names of interest (default=None)
            
        Adds
        ----------
        cell_scores : pd.DataFrame
            Cell count measurements
        """ 
        
        cts = np.copy(self.counts)
        names = self.names_cells
        
        df = pd.DataFrame(zip(names,
                              np.max(cts,axis=1),
                              np.sum(cts,axis=1)),
                          columns = ["names","max_count_value","total_counts"])
        df["counts_above0"] = np.count_nonzero(cts,axis=1)

        count1 = np.copy(cts)
        count1[np.where(count1 <= 1)] = 0 
        df["counts_above1"] = np.count_nonzero(count1,axis=1)
        
        count10 = np.copy(cts)
        count10[np.where(count10 <= 10)] = 0 
        df["counts_above10"] = np.count_nonzero(count10,axis=1)
        
        df["unique_count_values"] = pd.DataFrame(cts.T).nunique() - 1 # subtract 1 cause 0 doesnt count
        df["count_index"] = np.apply_along_axis(self._CountIndex,1,cts)
        
        if(gene_string is not None):
            ng = pd.Series(self.names_genes)
            gene_string_match = ng.str.match(gene_string)
            new_cts = cts[:,np.where(gene_string_match)[0]]
            df["feature_counts"] = np.sum(new_cts,axis=1)
        
        self.cell_scores = df
        
    def _SubsampleRow(self,row,n_counts):
        """
        Internal function for subsampling a row
        """ 
        
        transcripts = np.repeat(range(len(row)),row)
            
        new_row = np.zeros_like(row)
        unique,counts = np.unique(np.random.choice(transcripts,n_counts,replace=False),return_counts=True)
        np.put(new_row,unique,counts)
            
        return(new_row)
    
    def Subsample(self,n_counts):
        """
        Subsamples cells to a standard number of counts by randomly sampling observations without replacement.
        ...
        
        Parameters
        ----------
        n_counts : int
            The target number of counts, must be larger than the minimum total counts per cell
    
        
        Adds
        ----------
        subsample : numpy.ndarray
            New count matrix with row (cell) sums equal
        """ 
        
        depth_cells = np.sum(self.counts,axis=1)        
        assert all(np.sum(self.counts,axis=1) > n_counts), "subsample size must be larger than smallest count per cell"
        
        self.subsample = np.apply_along_axis(self._SubsampleRow,1,self.counts,n_counts=n_counts)

    def ScoreGenes(self,subsample=True):
        """
        Calculates several scores for count-based gene expression.
        ...
        
        Parameters
        ----------
        subsample : bool
            if true (default), score genes using subsampled counts
            otherwise, score counts using the unsampled count matrix
   
        Adds
        ----------
        gene_scores : pd.DataFrame
            Gene expression measurements
        """ 
        
        if (subsample == False):
            sg = self.counts
        else:
            assert hasattr(self,"subsample"), "expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix"
            sg = self.subsample
        
        names = self.names_genes
        
        df = pd.DataFrame(zip(names,
                              np.max(sg,axis=0),
                              np.sum(sg,axis=0)),
                          columns = ["names","max_count_value","total_counts"])
        df["counts_above0"] = np.count_nonzero(sg,axis=0)

        count1 = np.copy(sg)
        count1[np.where(count1 <= 1)] = 0 
        df["counts_above1"] = np.count_nonzero(count1,axis=0)
        
        count10 = np.copy(sg)
        count10[np.where(count10 <= 10)] = 0 
        df["counts_above10"] = np.count_nonzero(count10,axis=0)
        
        df["unique_count_values"] = pd.DataFrame(sg).nunique() - 1 # subtract 1 cause 0 doesnt count
        df["count_index"] = np.apply_along_axis(self._CountIndex,0,sg)
        
        self.gene_scores = df
        
    def Dot(self):
        """
        Calculates pairwise dot product between all cells based on counts   
        ...
        
        Adds
        ----------
        dots : numpy.ndarray
            similarity matrix of dot products between cells
        """         
        
        logging.info("Calculating dot products between rows...")
        n = self.counts.shape[0] # cells assumed to be rows

        #dots = self.counts @ self.counts.T
        dots = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                dot = np.dot(self.counts[i,:],self.counts[j,:])
                dots[i,j] = dot
                dots[j,i] = dot

        self.dots = dots
        
        logging.info("    done.")

    def Cluster(self,n_clusters):
        """
        Performs spectral clustering on dot products using sklearn  
        ...
        
        Parameters
        ----------
        n_clusters: int
            dimension of projection space
        
        Adds
        ----------
        cluster_labels : numpy.ndarray
            array of assigned cluster labels for each cell
        """ 
        
        logging.info("performing spectral clustering on dot products...")
        spectral = SpectralClustering(n_clusters=n_clusters, affinity='precomputed').fit(self.dots)
        self.cluster_labels = spectral.labels_
        logging.info("    done.")
        
    def PlotEmbedding(self):
        """
        Plot cells using spectral embedding of dot products. Clustering results and total counts are displayed
        """ 
        
        se = SpectralEmbedding(n_components=2, affinity='precomputed').fit(self.dots)
        self.spectral_embedding = se.embedding_
        
        fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(8,4))

        sns.scatterplot(ax = axes[0],
                        x = self.spectral_embedding[:,0],
                        y = self.spectral_embedding[:,1],
                        hue = self.cluster_labels, 
                        s = 10, 
                        palette = "tab10", 
                        linewidth=0) 
        axes[0].legend(loc=(1.04,0))
        axes[0].xaxis.set_ticklabels([])
        axes[0].yaxis.set_ticklabels([])
        axes[0].set(title="clusters")
        sns.scatterplot(ax = axes[1], 
                        x = self.spectral_embedding[:,0], 
                        y = self.spectral_embedding[:,1], 
                        hue = np.sum(self.counts,axis=1), 
                        s = 10, 
                        palette = "viridis", 
                        linewidth=0) 
        axes[1].legend(loc=(1.04,0))
        axes[1].xaxis.set_ticklabels([])
        axes[1].yaxis.set_ticklabels([])
        axes[1].set(title="total counts")

        fig.tight_layout()

    def RankMarkerGenes(self,method='prop-zero',subsample=True):
        """
        Ranks the top marker gene for each cluster from spectral clustering
        ...

        Parameters
        ----------
        method: str
            Rank by the difference in the proportion of cluster and non-cluster
            cells that have zero counts ('frac-zero' = default) or 
            the Wilcoxon rank sums statistic ('rank-sums')
        subsample : bool
            if true (default), score genes using subsampled counts
            otherwise, score counts using the unsampled count matrix
            
        Adds
        ----------
        marker_genes : pandas.DataFrame
            top 10 marker genes for each cluster
        marker_full : pandas.DataFrame
            list of dataframes of all ranked genes for each cluster 
        """ 
        if (subsample == False):
            sg = self.counts
        else:
            assert hasattr(self,"subsample"), "expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix"
            sg = self.subsample
            
        cluster = np.unique(self.cluster_labels)
        n = len(self.names_genes)
        
        rank_dfs = []
        top = pd.DataFrame()

        for c in cluster:
            cluster_cells = np.where(self.cluster_labels==c)
            noncluster_cells = np.where(self.cluster_labels!=c)

            if(method=="rank-sums"):
                res = np.zeros((n,2))
                for i in range(n):
                    res[i,:] = mannwhitneyu(self.counts[noncluster_cells,i].ravel(),self.counts[cluster_cells,i].ravel(),alternative="greater")
            
                df = pd.DataFrame(res,columns=['test-statistic','p-value']) 
                df['adjusted-p-value'] = multitest.fdrcorrection(res[:,1])[1]
                df['significant diff.'] = multitest.fdrcorrection(res[:,1])[0]
                df['rank'] = df['test-statistic'].rank(ascending=True).astype(np.uint)
            
            elif(method=="prop-zero"):
                cluster_positive = np.count_nonzero(self.counts[cluster_cells[0],:],axis=0)
                noncluster_positive = np.count_nonzero(self.counts[noncluster_cells[0],:],axis=0)
        
                prop_cluster = cluster_positive / np.size(cluster_cells)
                prop_noncluster = noncluster_positive / np.size(noncluster_cells)
                prop_a = prop_cluster - prop_noncluster
                
                df = pd.DataFrame(prop_a,columns=['diff. proportion zero']) 
                df['rank'] = df['diff. proportion zero'].rank(ascending=False).astype(np.uint)
            
            else:
                print ("error: unrecognized method")
                return
            
            df['gene index'] = range(len(self.names_genes))
            df['gene names'] = self.names_genes
            df['cluster_label'] = c
            rank_dfs.append(df)
            top = pd.concat((top,df.sort_values(by="rank").head(10)))
       
        self.marker_full = rank_dfs
        self.marker_genes = top
       
    def PlotMarker(self,gene_index):
        """
        Plot cells using spectral embedding and display counts in a given gene
        ...
        
        Parameters
        ----------
        gene_index: int
            index value of gene to be displayed        
        """ 
        marker_counts = self.counts[:,gene_index].ravel()

       # cmap = cm.get_cmap("inferno").copy()
       # cmap.set_under(color='#CFCFCF') 

        fig, axes = plt.subplots(1,2,figsize=(8,4))

        sns.scatterplot(ax = axes[0], 
            x = self.spectral_embedding[:,0],
            y = self.spectral_embedding[:,1], 
            hue = self.cluster_labels, 
            palette = "tab10", 
            s = 10,
            linewidth=0) 
        axes[0].legend(loc=(1.04,0))
        axes[0].xaxis.set_ticklabels([])
        axes[0].yaxis.set_ticklabels([])
        axes[0].set(title="clusters")
        sns.scatterplot(ax = axes[1],
            x = self.spectral_embedding[:,0], 
            y = self.spectral_embedding[:,1], 
            color = "#CFCFCF",
            s = 5,
            linewidth=0)
        sns.scatterplot(ax = axes[1],
            x = self.spectral_embedding[marker_counts > 0,0], 
            y = self.spectral_embedding[marker_counts > 0,1], 
            hue = marker_counts[marker_counts > 0],
            palette = "viridis",
            s = 10,
            linewidth=0) 
        axes[1].legend(loc=(1.04,0))
        axes[1].xaxis.set_ticklabels([])
        axes[1].yaxis.set_ticklabels([])
        axes[1].set(title="marker gene counts")
    
        fig.tight_layout()

    def PlotIMAElbow(self,max_features,u_bounds,subsample=True):
        """
        Plots the difference between the observe and reconstructed count matrix using integer matrix approximation and a series of total features
        ...
        
        Parameters
        ----------
        max_features: int
            maximum number of features to assess 
        u_bounds: list of int
            the upper bound for matrix U and matrix V
        subsample: bool
            if true (default), use the matrix of subsampled counts 
            otherwise, score counts using the unsampled count matrix
        
        """ 
        
        if (subsample == False):
            sg = self.counts
        else:
            assert hasattr(self,"subsample"), "expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix"
            sg = self.subsample
        
        norms = np.zeros(max_features)
        obs_norm = np.linalg.norm(sg,'fro')
        
        logging.info("approximating count matrix using up to %s features...",max_features)
        for i in range(max_features):
            params = IMA_params(i+1,u_bounds)
            U,V,Lambda_ = IMA(sg,params)
            norm = np.linalg.norm(U@np.diag(Lambda_)@V.T,'fro')
            norm_diff = obs_norm - norm
            norms[i] = norm_diff
        self.ima_differences = norms
        
        g = sns.scatterplot(x=range(1,max_features + 1),y = norms, color="black")
        g.set(xlabel="features",ylabel="diff. in matrix norms") 
        logging.info("    done.")

    def IMA(self,features,u_bounds,l_bounds=[0,0],maxiter=1000000.0,stop_crit=0.0001,subsample=True):
        """
        Performs integer matrix approximation on the count matrix
        ...
        
        Parameters
        ----------
        features: int
            target number of features 
        u_bounds: list of int
            the upper bound for matrix U and matrix V
        u_bounds: list of int
            the lower bound for matrix U and matrix V (default is 0 and 0)
        maxiter: float
            maximum number of iterations (default is 1000000)
        stop_crit: float
            criterion for stopping, based on difference between iterations 
        subsample: bool
            if true (default), use the matrix of subsampled counts 
            otherwise, score counts using the unsampled count matrix
            
        Adds
        ----------
        matrixU : numpy.ndarray
            matrix of dimensions cells x features
        matrixV : numpy.ndarray
            matrix of dimensions genes x features
        matrixLambda: numpy.ndarray
            diagonal scaling matrix
        """ 
        if (subsample == False):
            sg = self.counts
        else:
            assert hasattr(self,"subsample"), "expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix"
            sg = self.subsample

        logging.info("approximating count matrix...")
        params = IMA_params(features,u_bounds,l_bounds,maxiter,stop_crit)
        U,V,Lambda_ = IMA(sg,params)
    
        self.matrixU = U
        self.matrixV = V
        self.matrixLambda = np.diag(Lambda_)
        logging.info("    done.")

    def PlotIMA(self,x = 0, y = 1,subsample=True):
        """
        Plot cells using integer matrix approximation
        ...
        
        Parameters
        ----------
        x: int
           feature to plot on x-axis         
        y: int
           feature to plot on y-axis 
        subsample: bool
            if true (default), use the matrix of subsampled counts 
            otherwise, score counts using the unsampled count matrix
        
        """ 
        
        if (subsample == False):
            sg = self.counts
        else:
            assert hasattr(self,"subsample"), "expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix"
            sg = self.subsample
        
        loading = sg  @ (self.matrixV @ self.matrixLambda)

        fig, axes = plt.subplots(1,1,figsize=(4,4))

        sns.scatterplot(ax = axes,
            x = loading[:,x], 
            y = loading[:,y], 
            hue = self.cluster_labels,
            palette = "tab10", 
            s = 5,
            linewidth=0)
        axes.legend(loc=(1.04,0))
        axes.xaxis.set_ticklabels([])
        axes.yaxis.set_ticklabels([])
        axes.set(title="clusters")
        axes.set(xlabel = "".join(("feature: ",str(x))))
        axes.set(ylabel = "".join(("feature: ",str(y))))
        
        fig.tight_layout()
        
    def SharedCounts(self,n_clusters,n_cells=100,subsample=True):
        """
        Combines groups of genes with similar counts by clustering and summing
        ...
        
        Parameters
        ----------
        n_clusters: int
            number of clusters of genes 
        n_cells: int
            number of cells to sample for gene clustering
        subsample: bool
            if true (default), use the matrix of subsampled counts 
            otherwise, score counts using the unsampled count matrix
            
        Adds
        ----------
        sum_sharedcounts : numpy.ndarray
            similarity matrix of shared counts between genes
        sum_sharedcounts : numpy.ndarray
            matrix with counts summed within clusters of genes 
        sum_sharedcounts_all : numpy.ndarray
            same as sum_sharedcounts, but including genes not present in any cluster

        """ 
        
        if (subsample == False):
            sg = self.counts
        else:
            assert hasattr(self,"subsample"), "expecting array of subsampled counts, use subsample() or select subsample=False to use unsampled count matrix"
            sg = self.subsample
           
        logging.info("clustering genes based on shared counts...")
        # pull random cells and ensure they have at least one count    
        filtX = sg[np.random.randint(0,sg.shape[0],n_cells),:]
        filt_gene_index = np.where(np.sum(filtX,axis=0) > 0)[0] 
        filtX = filtX[:,filt_gene_index] 
        
        # calculate pairwise manhattan distances, convert to shared counts
        manhat = pairwise.manhattan_distances(filtX.T,sum_over_features = True)
        sumgenes = np.sum(filtX,axis = 0) 
        sumgenes_mat = np.add.outer(sumgenes,sumgenes) 
        self.sharedcounts = (sumgenes_mat.astype(np.uint) - manhat) / 2 
        
        # cluster based on shared counts
        spectral = SpectralClustering(n_clusters=n_clusters,affinity='precomputed').fit(self.sharedcounts) # cluster

        # group genes by clusters and sum counts
        combX = sg[:,filt_gene_index] 
        restX = np.delete(sg,filt_gene_index,1) # genes not in any cluster

        sumX = np.zeros((combX.shape[0],n_clusters),dtype=np.uint) 
        for s in range(n_clusters):
            sumX[:,s] = np.sum(combX[:,spectral.labels_ == s],axis=1)

        self.sum_sharedcounts = sumX[:,np.argsort(np.sum(sumX,axis=0))[::-1]]
        self.sum_sharedcounts_all = np.column_stack((self.sum_sharedcounts,restX))
        logging.info("    done.")
        
    def PlotSharedCounts(self,x = 0, y = 1):
        """
        Plot cells using matrix of counts summed by clusters of genes
        ...
        
        Parameters
        ----------
        x: int
           gene cluster to plot on x-axis         
        y: int
           gene cluster to plot on y-axis 

        """ 
            
        loading = self.sum_sharedcounts

        fig, axes = plt.subplots(1,1,figsize=(4,4))

        sns.scatterplot(ax = axes,
            x = loading[:,x], 
            y = loading[:,y], 
            hue = self.cluster_labels,
            palette = "tab10", 
            s = 5,
            linewidth=0)
        axes.legend(loc=(1.04,0))
        axes.xaxis.set_ticklabels([])
        axes.yaxis.set_ticklabels([])
        axes.set(title="clusters")
        axes.set(xlabel = "".join(("gene group: ",str(x))))
        axes.set(ylabel = "".join(("gene group: ",str(y))))
        
        fig.tight_layout()

    def _Normalize(self):
        """
        Recapitulates scanpy normalization
        ...
    
        Adds
        ----------
        norm_cells : numpy.ndarray
            normalization factor for each cell
        """        
        
        print("Normalizing counts...")
        # calculate normalization factors as in scanpy.pp.normalize_total
        self.norm_factor = 10000 / self.counts.sum(axis=1) 
        self.norm_counts = self.counts * self.norm_factor[:, np.newaxis]
        print("    done.")

    def _Log(self):
        """
        Recapitulates scanpy log transformation
        ...
    
        Adds
        ----------
        log_counts : list
            log transformed count matrix
        """        
        # logarithmize as in scanpy.pp.log1p
        print("Log transforming normalized counts...")
        self.log_counts = np.log( self.norm_counts + 1 )
        print("    done.")
        
    def _RescaleVariance(self):
        """
        Recapitulates scanpy scaling to unit variance
        ...
    
        Adds
        ----------
        scaled_counts : list
            counts scaled to gene unit variance
        """        
        # logarithmize as in scanpy.pp.log1p
        print("Scaling transformed counts to gene unit variance...")
        self.scaled_counts = preprocessing.scale(self.log_counts,axis=0,with_mean = False)
        print("    done.")  
        
    def _Center(self):
        """
        Recapitulates scanpy centering scaled and transformed data
        ...
    
        Adds
        ----------
        centered_counts : list
            centers counts to zero
        """        
        # logarithmize as in scanpy.pp.log1p
        print("Centering counts...")
        self.centered_counts = preprocessing.scale(self.scaled_counts,axis=0,with_std = False)
        print("    done.")  