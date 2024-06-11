import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, davies_bouldin_score, silhouette_score
import time
from datetime import datetime
import os
import h5py
from .preprocess import Scaler

class vqpca:

    """
    Vector Quantization Principal Component Analysis (VQPCA) is a clustering 
    algorithm that combines vector quantization with principal component analysis.
    It partitions data into clusters such that the local reconstruction
    error is minimized. This class provides methods for fitting the VQPCA 
    model to data, reconstructing data, clustering new data, evaluating model 
    performance, and writing output files.

    Attributes:
        X_ : array-like
            Data matrix to be analyzed.
        n_obs_ : int
            Number of observations in the data matrix.
        nvars_ : int
            Number of variables in the data matrix.
        k_ : int
            Number of clusters.
        q_ : float
            Variance threshold for stopping rule in PCA.
        stopping_rule_ : str
            Rule to determine when to stop the PCA algorithm.
        itmax_ : int
            Maximum number of iterations for convergence.
        atol_ : float
            Absolute tolerance for convergence.
        rtol_ : float
            Relative tolerance for reconstruction error.
        ctol_ : float
            Tolerance for centroids.
        pca_ : list of PCA objects
            PCA objects for each cluster.
        labels_ : array-like
            Labels indicating the cluster to which each point belongs.
        eps_ : float
            Reconstruction error.
        C_ : array-like
            Centroids of the clusters.

    Methods:
        __init__() : Initialize the vqpca object.
        initialize_centroids() : Initialize the centroids of the clusters.
        initialize_pca() : Initialize PCA objects for each cluster.
        fit() : Fit the VQPCA model to the data.
        reconstruct() : Reconstruct the data using the fitted VQPCA model.
        predict() : Cluster a new matrix based on the local basis found by local PCA.
        get_components() : Get the principal components from each local PCA.
        write_output_files() : Create output files for the VQPCA model.
        evaluate() : Evaluate the VQPCA model using a specified evaluation metric.
    """

    def __init__(self, X, k=2, stopping_rule="variance", q=0.99, itmax=200, 
                 atol=1e-8, rtol=1e-8, ctol=1e-6):
        
        """
        Initialize the VQPCA with data and parameters.

        Parameters
        ----------
        X : array-like
            The data matrix to be analyzed.
        k : int, optional
            The number of clusters (default is 2).
        stopping_rule : str, optional
            The rule to determine when to stop the PCA algorithm (default is "variance").
        q : float, optional
            The variance threshold for the stopping rule (default is 0.99).
        itmax : int, optional
            The maximum number of iterations for convergence (default is 200).
        atol : float, optional
            The absolute tolerance for convergence (default is 1e-8).
        rtol : float, optional
            The relative tolerance for reconstruction error (default is 1e-8).
        ctol : float, optional
            The tolerance for centroids (default is 1e-6).

        Raises
        ------
        ValueError
            If the specified stopping rule does not exist.
        """

        # Data matrix
        self.X_ = X
        (nobs, nvars) = np.shape(X)
        self.n_obs_ = nobs
        self.nvars_ = nvars

        # Number of clusters
        self.k_ = k

        # PCA stopping input
        self.q_ = q

        # Stopping rule (check existence first)
        stop_list = ["variance", "n_eigs"]
        if stopping_rule not in stop_list:
            raise ValueError("The stopping rule specified does not exist")
        self.stopping_rule_ = stopping_rule

        # Parameters for convergence
        self.itmax_ = itmax
        self.atol_  = atol  
        self.rtol_  = rtol  # Tolerance on reconstruction error
        self.ctol_  = ctol  # Tolerance on centroids


    # Initialization methods
    def initialize_centroids(self, method, Ci=None):

        """
        Initialize the cluster centroids using a specified method.

        Parameters
        ----------
        method : str
            The method to initialize the centroids. Available methods are: 
            'uniform', 'kmeans', 'random'.
        Ci : array-like, optional
            Initial centroids, used only if method supports external input 
            (currently not used in this implementation).

        Returns
        -------
        self : vqpca
            The instance of the class with initialized centroids.

        Raises
        ------
        ValueError
            If the specified initialization method is not recognized.
        """
                
        if method == 'kmeans':
            # Centroids are initialized using k-means
            kmeans = KMeans(n_clusters=self.k_, random_state=0).fit(self.X_)
            self.labels_ = kmeans.labels_
            self.C_ = kmeans.cluster_centers_

        elif method == 'random':
            np.random.seed(1000)
            # Centroids are initialized randomly
            self.labels_ = np.random.randint(0, self.k_, self.X_.shape[0])
            self.C_ = np.zeros((self.k_, self.X_.shape[1]))
            for i in range(self.k_):
                self.C_[i,:] = np.mean(self.X_[self.labels_==i,:], axis=0)

        elif method == 'uniform':
            # Centroids are initialized uniformly
            npts = self.X_.shape[0]
            # Get number of groups
            ngroups = int(np.floor(npts / self.k_))
            # Get number of points in the last group
            nlast = npts - ngroups * self.k_
            # Initialize labels
            self.labels = np.zeros(npts)
            # Initialize centroids
            self.C_ = np.zeros((self.k_, self.X_.shape[1]))
            # Loop over groups
            for i in range(self.k_):
                # Assign labels
                self.labels[i*ngroups:(i+1)*ngroups] = i
                # Assign centroids
                self.C_[i,:] = np.mean(self.X_[i*ngroups:(i+1)*ngroups,:], axis=0)
            # Assign labels to last group
            self.labels[-nlast:] = self.k_ - 1
            # Assign centroids to last group
            self.C_[-1,:] = np.mean(self.X_[-nlast:,:], axis=0)

        else:
            raise ValueError("Initialization method for centroids not recognized")

        return self
        
    def initialize_pca(self, n_start=2):

        """
        Initialize PCA objects for each cluster.

        Parameters
        ----------
        n_start : int, optional
            The number of principal components to start with for each PCA 
            (default is 2).

        Returns
        -------
        self : vqpca
            The instance of the class with initialized PCA objects.
        """

        pca_list = [] 
        for i in range(self.k_):
            pca = PCA(n_components=n_start)
            idps = np.where(self.labels_==i)[0]
            pca.fit(self.X_[idps,:])
            pca_list.append(pca)

        self.pca_ = pca_list
        return self

    def fit(self, n_start=2, verbose=True, init_centroids='kmeans'):

        """
        Fit the VQPCA model to the data.

        Parameters
        ----------
        n_start : int, optional
            The number of principal components to start with for each PCA 
            (default is 2).
        verbose : bool, optional
            Whether to print progress information (default is True).
        init_centroids : str, optional
            The method to initialize centroids ('kmeans', 'random', 'uniform').
            (default is 'kmeans').

        Returns
        -------
        labels_new : array-like
            The labels of the data points after clustering.
        self : vqpca
            The instance of the class with fitted model.

        Notes
        -----
        This method fits the VQPCA model to the data using the specified parameters.

        """ 

        # Initialize global time
        st_global = time.time()

        # Initialize flag for convergence
        conv = False
        # Initialize iteration counter
        iter = 0

        # Initialize centroids
        self.initialize_centroids(method=init_centroids)
        
        # Initialize pca
        self.initialize_pca(n_start=n_start)

        while conv == False and iter < self.itmax_:

            # Monitor iteration time
            st_iter = time.time()
            # Matrix of projection errors
            eps_matrix = np.zeros((self.X_.shape[0], self.k_))
            # Perform PCA in each cluster and calculate projection error
            for i in range(self.k_):
                # Predict the point using the local PLS model
                U_scores = self.pca_[i].transform(self.X_)
                # Reconstruct X
                X_rec    = self.pca_[i].inverse_transform(U_scores)
                # Calculate norm of the residual
                eps_matrix[:,i] = np.sum((self.X_ - X_rec)**2, axis=1)
            # Assign each point to the cluster with the lowest residual
            labels_new = np.argmin(eps_matrix, axis=1)
            # Calculate the new centroids
            centroids_new = np.zeros((self.k_, self.X_.shape[1]))
            for i in range(self.k_):
                centroids_new[i,:] = np.mean(self.X_[labels_new==i,:], axis=0)
            # Calculate the change in the centroids
            delta_centroids = np.linalg.norm(centroids_new - self.C_) / np.linalg.norm(self.C_)
            if delta_centroids < self.ctol_:
                conv = True
                if verbose:
                    print('Converged because of the centroids')

            # Check convergence on reconstruction error variance
            eps_rec_new = np.min(eps_matrix, axis=1)
            if iter == 0:
                self.eps_ = np.mean(eps_rec_new)
            else:
                delta_eps = np.abs((self.eps_ - np.mean(eps_rec_new)) / self.eps_)
                if delta_eps < self.rtol_:
                    conv = True
                    if verbose:
                        print('Converged because of reconstruction error variance')

            # Update PCA objects
            pca_list = []
            rec_err_clust = []
            for i in range(self.k_):
                pca = PCA(n_components=self.q_)
                pca.fit(self.X_[labels_new==i,:])
                pca_list.append(pca)
                U_local = pca.transform(self.X_[labels_new==i])
                X_local_rec = pca.inverse_transform(U_local)
                rec_err_clust.append(mean_squared_error(self.X_[labels_new==i], X_local_rec))
                
            # Update self
            self.pca_ = pca_list
            self.C_ = centroids_new
            self.labels_ = labels_new
            self.eps_    = np.mean(eps_rec_new)

            # Update iteration counter 
            iter += 1
            if iter == self.itmax_:
                print("Iterations reached maximum allowed number")

            # Time for iteration
            et_iter = time.time()
            dt_iter = et_iter - st_iter

            if verbose == True:
                print("Iteration: ", iter)
                print("Mean squared global reconstruction error = ", self.eps_ )
                print("Mean squared error in each cluster: ")
                print(rec_err_clust)
                if iter > 1:
                    print("Reconstruction error variance = ", delta_eps)
                    print("Centroids variance = ", delta_centroids)
                print(f"Elapsed time for iteration: {dt_iter} seconds \n \n")

        # Measure time at convergence
        et_global = time.time()
        dt_global = et_global - st_global
        if verbose == True:
            print(f"Elapsed time for global convergence: {dt_global} seconds")

        # Display information
        if verbose:
            print("Obtained reconstruction error = ", self.eps_)

        return labels_new, self
    
    def reconstruct(self, X=None):

        """
        Reconstruct the data using the fitted VQPCA model.

        Parameters
        ----------
        X : array-like, optional
            The data to be reconstructed. If not provided, the original data 
            used for fitting will be reconstructed (default is None).

        Returns
        -------
        X_rec : array-like
            The reconstructed data.

        Notes
        -----
        This method provides the reconstructed data using the fitted VQPCA model. 
        If the input data X is not provided, it reconstructs the original data 
        used for fitting.

        """

        # Check if X was given
        if X.all() == None:
            X = self.X_

        # Initialize array of reconstructed X
        X_rec = np.zeros_like(X)
        for i in range(self.k_):
            # Get local PCA scores
            U_local = self.pca_[i].transform(X[self.labels_==i])
            # Reconstruct back using inverse transformation        
            X_rec[self.labels_==i] = self.pca_[i].inverse_transform(U_local)

        return X_rec
    
    def predict(self, X):

        """
        Cluster a new matrix X based on the local basis found by local PCA.

        Parameters
        ----------
        X : array-like
            The new matrix to be clustered.

        Returns
        -------
        labels : array-like
            The labels indicating the cluster to which each point in X belongs.

        Raises
        ------
        ValueError
            If the number of columns in X is different from the original data.

        Notes
        -----
        This method clusters a new matrix X based on the local basis that 
        the local PCA found in each cluster.

        """
        
        # Check shape of X
        (nobs, nvars) = np.shape(X)
        if nvars != self.nvars_:
            raise ValueError("Number of columns of new matrix must be the same as original")
        
        # Matrix of projection errors
        eps_matrix = np.zeros((X.shape[0], self.k_))
        # Perform PCA in each cluster and calculate projection error
        for i in range(self.k_):
            # Predict the point using the local PLS model
            U_scores = self.pca_[i].transform(X)
            # Reconstruct X
            X_rec    = self.pca_[i].inverse_transform(U_scores)
            # Calculate norm of the residual
            eps_matrix[:,i] = np.sum((X - X_rec)**2, axis=1)

        # Assign each point to the cluster with the lowest residual
        labels = np.argmin(eps_matrix, axis=1)

        return labels

    ### Function for getting access to attributes ###
    def get_components(self):

        """
        Get the principal components from each local PCA.

        Returns
        -------
        components : list of arrays
            The principal components from each local PCA.

        Notes
        -----
        This function returns the principal components obtained from each local PCA model.
        Each element in the list corresponds to the principal components of a cluster.

        """

        components = []
        for i in range(self.k_):
            ci = self.pca_[i].components_.T
            components.append(ci)
        
        return components
    
    def write_output_files(self, foldername=None, format='h5'):

        """
        Create output files for the VQPCA model.

        Parameters
        ----------
        foldername : str, optional
            The name of the folder to save the output files. If None, a folder 
            with a timestamp will be created (default is None).
        format : str, optional
            The format of the output files ('h5' or 'txt'). 'h5' is recommended 
            for large datasets, 'txt' is suitable for smaller datasets 
            (default is 'h5').

        Returns
        -------
        self : vqpca
            The instance of the class.

        Notes
        -----
        This method creates output files containing components, labels, 
        reconstruction error, and centroids of the VQPCA model.

        """
        
        if foldername == None:
            # Get the current date
            current_date = datetime.now().strftime("%Y-%m-%d")
            foldername = "VQPCA_k" + str(self.k_) + "q_" + str(self.q_) + '_' + current_date
        
        # Create directory if it does not exists
        if os.path.exists(foldername) == False:
            os.mkdir(foldername)

        if format == "h5":
            file_path = foldername + '/output.h5'
            with h5py.File(file_path, 'w') as f:

                # Write components
                components = self.get_components()
                basis_group = f.create_group("components")
                for i, array in enumerate(components):
                    basis_group.create_dataset(f'array_{i}', data=array)

                # Write labels
                f.create_dataset("labels", data=self.labels_)

                # Reconstruction error
                f.create_dataset("reconstruction_error", data=self.eps_)

                # Write centroids
                f.create_dataset("centroids", data=self.C_)

        elif format == "txt":
            
            np.savetxt(foldername + '/labels.txt', self.labels_, fmt='%d')
            np.savetxt(foldername + '/centroids.txt', self.C_)

        return self

    ### VQPCA evaluation ###
    def evaluate(self, score='ILPCA'):
        """
        Evaluate the VQPCA model using a specified evaluation metric.

        Parameters
        ----------
        score : str, optional
            The evaluation metric to use ('ILPCA', 'DB', 'silhouette').
            (default is 'ILPCA').

        Returns
        -------
        metric : float
            The evaluation metric score.

        Notes
        -----
        This method evaluates the VQPCA model using the specified evaluation metric.

        """

        # ILPCA score evaluation
        if score == 'ILPCA':

            # Calculate squared reconstruction error in the clusters
            rec_err_clust = np.zeros(self.k_)
            for i in range(self.k_):
                U_scores = self.pca_[i].transform(self.X_[self.labels_==i])
                X_rec = self.pca_[i].inverse_transform(U_scores)
                rec_err_clust[i] = (np.mean(np.sum((self.X_[self.labels_==i] - X_rec)**2, axis=1)))

            # Initialize metric
            metric = 0.0
            for i in range(self.k_):
                # Reconstruction error in cluster i
                eps_i = rec_err_clust[i]
                # Initialize db
                db_iter = 0.0
                for j in range(self.k_):
                    if j != i:
                        # Reconstruction error in cluster j
                        eps_j = rec_err_clust[j]
                        # Merge cluster i and j
                        X_ij = np.vstack((self.X_[self.labels_==i,:], self.X_[self.labels_==j,:]))
                        # Perform PCA in the merged cluster
                        pca = PCA(n_components=self.q_)
                        pca.fit(X_ij)
                        # Reconstruct Xij
                        U_scores = pca.transform(X_ij)
                        X_rec    = pca.inverse_transform(U_scores)
                        # Reconstruction error of merged cluster
                        eps_ij = np.mean(np.sum((X_ij - X_rec)**2, axis=1))
                        # Get max between all the clusters pairs
                        db_iter = max(db_iter, (eps_i + eps_j)/eps_ij)

                metric += db_iter # Sum for the clusters
            
            # Average
            metric = metric / self.k_

        # If DB index was chosen
        elif score == "DB":
            metric = davies_bouldin_score(self.X_, self.labels_)

        elif score == "silhouette":
            metric = silhouette_score(self.X_, self.labels_)

        return metric

