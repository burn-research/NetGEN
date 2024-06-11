import numpy as np

# Class for scaler function similarly to sklearn
class Scaler:
    """
    A class used to scale data using various methods.

    Attributes
    ----------
    method : str
        The scaling method to be used ('auto', 'pareto', 'range', 'vast', 'max').
    center : bool, optional
        Whether to center the data before scaling (default is True).
    """

    # Default constructor
    def __init__(self, method, center=True):

        """
        Initialize the Scaler with a specified method and centering option.

        Parameters
        ----------
        method : str
            The scaling method to be used.
        center : bool, optional
            Whether to center the data before scaling (default is True).
        """

        self.method = method
        self.center = center

    # Fit function
    def fit(self, X):
        """
        Compute the mean (if centering) and scaling factors for the data.

        Parameters
        ----------
        X : array-like
            The data to compute the scaling factors for.

        Returns
        -------
        self : Scaler
            The fitted Scaler object.
        """
    
        # Number of columns
        ncols = np.shape(X)[1]
        # Check if center is true
        if self.center == True:
            c = np.mean(X,0)
        else:
            c = np.zeros(ncols)

        # Calculate the scaling value according to X
        if self.method == 'auto':
            gamma = np.std(X,0)
        elif self.method == 'pareto':
            gamma = np.sqrt(np.std(X,0))
        elif self.method == 'range':
            gamma = np.max(X,0)-np.min(X,0)
        elif self.method == 'vast':
            gamma = np.std(X,0)*np.mean(X,0)
        elif self.method == 'max':
            gamma = np.max(X,0)
        else:
            print('No method specified or the method was not recognized')
            gamma = np.ones(ncols)
        
        # Check if any of gamma is lower than 1e-16
        for i in range(len(gamma)):
            if abs(gamma[i]) < 1e-16:
                gamma[i] = 1e-16
                Warning('Division by zero avoided')

        # Assign c to the class
        self.c = c
        self.gamma = gamma

        return self
    
    # Fit_transform function
    def fit_transform(self, X, center=True):
        """
        Fit the Scaler to the data and then transform it.

        Parameters
        ----------
        X : array-like
            The data to be scaled.
        center : bool, optional
            Whether to center the data before scaling (default is True).

        Returns
        -------
        Xs : array-like
            The scaled data.
        """

        # Fit
        self.fit(X)
        # Transform
        if center == True:
            Xs = (X - self.c) / self.gamma
        else:
            Xs = X / self.gamma
        return Xs
    
    # Transform function
    def transform(self, X, center=True):

        """
        Scale the data using the previously computed mean and scaling factors.

        Parameters
        ----------
        X : array-like
            The data to be scaled.
        center : bool, optional
            Whether to center the data before scaling (default is True).

        Returns
        -------
        Xs : array-like
            The scaled data.

        Raises
        ------
        ValueError
            If the Scaler has not been fitted yet.
        """

        # Check if c and gamma exists already
        if hasattr(self, 'gamma') and hasattr(self, 'c') == False:
            raise ValueError("Fit the scaler before transforming other data using fit()")
        else:
            if center == True:
                Xs = (X - self.c) / self.gamma
            else:   
                Xs = X / self.gamma
        return Xs
    
    # Inverse transformation
    def inverse_transform(self, X):

        """
        Reverse the scaling of the data using the previously computed mean and scaling factors.

        Parameters
        ----------
        X : array-like
            The scaled data to be reversed.

        Returns
        -------
        Xi : array-like
            The original data before scaling.

        Raises
        ------
        ValueError
            If the Scaler has not been fitted yet.
        """

        # Check if c and gamma exists already
        if hasattr(self, 'gamma') and hasattr(self, 'c') == False:
            raise ValueError("Fit the scaler before transforming other data using fit()")
        Xi = (X + self.c) * self.gamma
        return Xi
    
