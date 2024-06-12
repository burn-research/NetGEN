import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

from .clustering import *
from .dataimport import *
from .graphs import *
from .preprocess import *
from .netsmoke import *

# Create class for CRN generation
class CRNGen:

    """
    CRNGen is a class for the automatic generation of chemical reactor networks
    from CFD data. It involves several steps:
    - Unsupervised Learning algorithms for clustering
    - Graph-based reassignment of data points
    - Computing mass flowrates and getting reactors' attributes (volume, temperature...)
    - Simulating CRNs (available interface with NetSMOKEpp solver)

    Parameters:
    -----------
    datainterface (str): Data interface type, default is 'fluent'.
    crninterface (str): CRN interface type, default is 'netsmoke'.
    """

    def __init__(self, datainterface='fluent', crninternface='netsmoke'):

        # Default constructur
        self.datainterface_ = datainterface
        self.crninterface_  = crninternface

    # Set datadictionary
    def setDataDictionary(self, datainfo, datapath, dim):

        """
        Set the data dictionary for the object.
        
        Parameters:
        -----------
        datainfo (dict): Dictionary containing data information.
        datapath (str): Path to the data.
        dim (int): Dimensionality of the data.
        
        Returns:
        --------
        self: The object instance with updated data dictionary.
        """

        # Check fields in datadictionary
        if 'volume' not in datainfo:
            raise Warning("The 'volume' field is not in datainfo")
        if 'density' not in datainfo:
            raise Warning("The 'density' field is not in datainfo")
        if 'solution' not in datainfo:
            raise Warning("The 'solution' field is not in datainfo")
        if 'velocity' not in datainfo:
            raise Warning("The 'velocity' field is not in datainfo")
        if 'angle' not in datainfo:
            raise Warning("The 'angle' field is not in datainfo")

        # Check what interface should be used
        if self.datainterface_ == 'fluent':
            datadictionary = GenDataFluent(datainfo, datapath, dim, verbose=False)
        else:
            raise ValueError("datainterface not recognized or not yet available")
        
        # Set attribute
        self.datadictionary_ = datadictionary
        self.dim_ = dim
        return self
    
    # Graph setter
    def setGraph(self, filepath=None, verbose=False, option='fluent',
                      connectivities=None, massflowrates=None):
        
        """
        Set the graph for the crngen object based on the specified option, namely
        the connectivities between mesh cells.
        
        Parameters:
        -----------
        filepath (str, optional): Path to the file containing graph data (required if option is 'fluent').
        verbose (bool): Flag to enable verbose output.
        option (str): Option for graph generation ('fluent' or 'custom').
        connectivities (numpy.ndarray, optional): Connectivity matrix (required if option is 'custom').
        massflowrates (numpy.ndarray, optional): Mass flow rates (optional for 'custom' option).
        
        Returns:
        --------
        self: The object instance with updated graph.
        """
        
        # Check for valid option
        if option != 'fluent' and option != 'custom':
            raise ValueError("Specified option for Graph generation was not recognized")
        
        # Check for filepath
        if option == 'fluent':
            if filepath == None:
                raise ValueError("Option fluent for graph generation was specified but no filepath was given")
        
        # Check for correct inputs for the custom method
        if option == 'custom':
            if connectivities is None:
                raise ValueError("Custom option was specified but no connectivty matrix was given")
            elif connectivities is not None and massflowrates is None:
                raise Warning("Custom option was specified but the mass flowrates were not given")
            
        # Set graph
        if option == 'fluent':
            G = GenGraph(filepath=filepath, option=option)
            # Import connectivity file
            dt  = pd.read_csv(filepath, sep='\t')
            v1  = dt.values[:,3].astype(int)
            v2  = dt.values[:,4].astype(int)
            mfi = dt.values[:,5]
            # Remove boundary cells
            mfi = mfi[v2 != -1].reshape(-1,1)
            v1  = v1[v2 != -1].reshape(-1,1)
            v2 = v2[v2 != -1].reshape(-1,1)
        
        elif option == 'custom':
            G = GenGraph(connectivities=connectivities, massflowrates=massflowrates)
            v1 = connectivities[:,0]
            v2 = connectivities[:,1]
            mfi = massflowrates

        # Set attributes
        self.G_ = G
        self.v1_ = v1
        self.v2_ = v2
        self.mfi_ = mfi

        return self
    
    def setData(self, dataset='state-space', X=None):

        """
        Set the data for the object based on the specified dataset type.
        
        Parameters:
        -----------
        dataset (str): Type of dataset to set ('state-space', 'reduced', or custom).
        X (numpy.ndarray, optional): Custom data to set (required if dataset is custom).
        
        Returns:
        --------
        self: The object instance with updated data.
        """

        if dataset == 'state-space':
            if 'solution' not in self.datadictionary_:
                raise ValueError("state-space was specified as setter for data but the 'solution' field is not self.datadictionary_")
            data = self.datadictionary_['solution']['values']

        elif dataset == 'reduced':
            # Extract temperature
            T = self.datadictionary_['solution']['values'][:,0].reshape(-1,1)
            # Extract velocity
            v = self.datadictionary_['velocity']['values'][:,0].reshape(-1,1)
            # Get residence time
            t = self.datadictionary_['tau']['values'][:,0].reshape(-1,1)
            # Process velocity
            v = v ** 1.5
            # Preprocess tau
            t = np.log10(t); t[t>1] = 1
            # Create data matrix
            data = np.concatenate((T, v, t), axis=1)

        else:
            if X is None:
                raise ValueError("Not a valid dataset method and X was not specified.")
            # Check shape
            if np.shape(X)[0] != len(self.datadictionary_['cell_ids']):
                raise ValueError("Number of rows of X is different from the number of cells in self.datadictionary_")    
            # If X is correctly given
            data = X

        # Set data as attribute
        self.data_ = data
        return self
    
    def RunClustering(self, nk, 
                      alg='k-means', graph_reassign=True, 
                      scaling='auto', center=True,
                      random_state=42, 
                      # VQPCA options
                      stop_rule="variance", q=0.99, cluster_itmax=300, init='random',
                      reassign_method='volume', v_thresh=1e-6, 
                      verbose=False,):
        
        """
        Run a clustering algorithm on the data.

        Parameters:
        -----------
        nk (int): Number of clusters.
        alg (str): Algorithm to use ('k-means' or 'vqpca').
        graph_reassign (bool): Flag to enable graph reassignment.
        scaling (str): Scaling method for the data.
        center (bool): Flag to enable centering of the data.
        random_state (int): Seed for random number generation.
        stop_rule (str): Stopping rule for vqpca.
        q (float): Parameter for vqpca.
        cluster_itmax (int): Maximum iterations for clustering.
        init (str): Initialization method for clustering.
        reassign_method (str): Method for graph reassignment.
        v_thresh (float): Threshold for volume reassignment.
        verbose (bool): Flag to enable verbose output.

        Returns:
        --------
        self: The object instance with updated clustering labels.
        """

        # Check if data is attribute
        if hasattr(self, 'data_') == False:
            raise ValueError("Set the data using the self.setData method first")
        elif graph_reassign == True and hasattr(self, 'G_') == False:
            raise ValueError("Set the Graph using the self.setGraph method first")
        
        # Check for option alg
        if alg != 'k-means' and alg != 'vqpca':
            raise ValueError("Specified unsupervised algorithm does not exist")
        
        # Scale the data
        scaler = Scaler(method=scaling, center=center)
        X_scaled = scaler.fit_transform(self.data_)

        # Run the unsupervised clustering algorithm
        if alg == 'k-means':
            kmeans = KMeans(n_clusters=nk, init=init, random_state=random_state, max_iter=cluster_itmax)
            kmeans.fit(X_scaled)
            labels = kmeans.labels_
        
        elif alg == 'vqpca':
            model = vqpca(X_scaled, k=nk, stopping_rule=stop_rule, itmax=cluster_itmax, q=q)
            labels = model.fit(verbose=verbose, init_centroids=init)
        
        # Now graph reassignment
        if graph_reassign:
            labels_new = ReassignNodes(self.G_, labels, self.datadictionary_, method=reassign_method, 
                                       v_thresh=v_thresh, verbose=verbose)
            
            # Match with number of reactors
            labels_new = MatchNClusters(self.G_, labels_new, self.datadictionary_, nk, verbose=verbose)
        else:
            labels_new = labels

        # Set attributes
        self.labels_ = labels_new
        self.nr_ = max(labels_new)+1

        return self
    
    def plotClusters(self, savefig=False, figname=None, figsize=(6,4)):

        """
        Plot the clustering solution on a scatter plot.

        Parameters:
        -----------
        savefig (bool): Flag to save the figure.
        figname (str, optional): Filename to save the figure.
        figsize (tuple): Size of the figure.

        Returns:
        --------
        None
        """

        xx = self.datadictionary_['coordinates']['values'][:,0]
        yy = self.datadictionary_['coordinates']['values'][:,1]

        fig, ax = plt.subplots(figsize=figsize)
        sc = ax.scatter(xx, yy, c=self.labels_, s=1)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_title('Clustering solution')
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label('Cluster index')
        fig.tight_layout()
        plt.show()
        if savefig:
            plt.savefig(figname, dpi=600)

    def getMassFlowrates(self):

        """
        Calculate and set the mass flow rates between clusters.

        Returns:
        self: The object instance with updated mass flow rate matrix.
        """

        # Get mass flowrates
        Mf = np.zeros((self.nr_, self.nr_))
        # Get left cells
        v1 = self.v1_
        v2 = self.v2_
        mfi = self.mfi_
        # Loop over connection to check neighbor clusters
        for i in range(len(v1)):
            # Left cell's id (coherent with Python)
            id1 = v1[i] - 1
            # Right cell's id (coherent with Python)
            id2 = v2[i] - 1
            # Elementary mass flowrate
            mfii = mfi[i]
            # Left cluster id
            c1 = self.labels_[id1]
            # Right cluster id
            c2 = self.labels_[id2]
            # If different update mass flowrate matrix
            if c1 != c2:
                # If positive Mf is going from c1 to c2
                if mfii > 0:
                    Mf[c1,c2] = mfii + Mf[c1,c2]
                # If negative is the opposite
                else:
                    Mf[c2,c1] =  - mfii + Mf[c2,c1]
        # Set attribute
        self.Mf_ = Mf
        return self

    def getInletStreams(self, streams, filepath=None, verbose=False):

        """
        Calculate and set the inlet and outlet mass flow rates based on the given streams.

        Parameters:
        streams (list of dict): List of streams with 'id', 'fluentid', 'name', and 'gas'.
        filepath (str, optional): Path to the boundary cell file (required if datainterface is 'fluent').
        verbose (bool): Flag to enable verbose output.

        Returns:
        self: The object instance with updated inlet and outlet mass flow rates.
        """

        if self.datainterface_ == 'fluent':
            if filepath is None:
                raise ValueError("Fluent data interface in use, you must specify the path to the boundary cell file. See Fluent UDF")
            
            # Import the file
            dt = pd.read_csv(filepath, sep='\t', header=None)
            bcs = dt.values

            # Initialize inlet mass flowrates
            inlets_mf = np.zeros(self.nr_)
            # Initialize inlet stream
            inlet_streams = np.zeros((self.nr_, len(streams)))
            # Initialize outlet mass fllowrates
            outlets_mf = np.zeros(self.nr_)

            # Scan through the boundary cells
            for i in range(np.shape(bcs)[0]):
                # Id of the cell coherent with the graph
                cell_id = int(bcs[i,0])
                # Id of the Fluent zone
                zone_id = bcs[i,1]
                # Mass flowrate across the cell
                mi      = bcs[i,2]
                # Get the cluster of the cell
                ci = self.labels_[cell_id-1]
                # If mi is positive, then the cell it's an outlet
                if mi > 0:
                    outlets_mf[ci] += mi
                # If negative, the cell it's an inlet
                else:
                    inlets_mf[ci] -= mi
                    # Check the respective zone id
                    for j in range(len(streams)):
                        # Get id of the stream (j-th column of array inlet streams)
                        stream_id = streams[j]['id']
                        # Check fluent zone of the stream 
                        stream_zone_id = streams[j]['fluentid']
                        # If zone_id and stream_zone_id are the same get mass flowrate
                        if stream_zone_id == zone_id:
                            inlet_streams[ci, stream_id] -= mi
            
            # Now check how many inlets have been detected
            if verbose:
                for i in range(self.nr_):
                    for j in range(len(streams)):
                        if inlet_streams[i,j] != 0:
                            print("Reactor ", i, " receive inlet from stream ", streams[j]['id'], " called ", streams[j]['name'])

        # Initialize list of inlet mixtures
        inlet_mixtures = [None] * self.nr_

        # Check for multiple inlets
        for i in range(self.nr_):
            count = 0
            idms  = []
            # Scan through inlets
            if inlets_mf[i]>0:
                for j in range(inlet_streams.shape[1]):
                    if inlet_streams[i,j]>0:
                        count += 1
                        idms.append(j)
                # If multiple inlet we mix the streams
                if count > 1:
                    print("Multiple inlet detected")
                    mxs = []
                    # Get gas objetcs
                    for j in range(len(idms)):
                        g = streams[j]['gas']   # Get gas object
                        q = ct.Quantity(g, mass=inlet_streams[i,idms[j]], constant="HP")
                        mxs.append(q)
                    M = mxs[0]
                    for j in range(len(mxs)-1):
                        M = M + mxs[j+1]
                    # Append to inlet mixtures
                    inlet_mixtures[i] = M

                # Single inlet
                else:
                    inlet_mixtures[i] = streams[idms[0]]['gas']

        # Update attributes
        self.inlet_mixtures_ = inlet_mixtures
        self.inlet_streams_ = inlet_streams
        self.outlets_mf_ = outlets_mf
        self.inlets_mf_  = inlets_mf
        self.streams_    = streams

        return self
    
    def checkMassBalance(self, tol=1e-2, verbose=False):

        """
        Check the mass balance of the system.

        Parameters:
        -----------
        tol (float): Tolerance level for imbalance.
        verbose (bool): Flag to enable verbose output.

        Returns:
        --------
        self: The object instance with updated mass balance attributes.
        """

        # Get inlet mass flowrates
        m_in  = self.inlets_mf_
        m_out = self.outlets_mf_

        # Get internal mass flowrates
        mf = self.Mf_

        # Check global balance
        geps = (np.sum(m_in) - np.sum(m_out))/max(np.sum(m_in), np.sum(m_out))
        if geps > tol:
            raise ValueError(f"Global imbalance of {geps} detected")
        
        # Check local imbalance
        leps = np.zeros(self.nr_)
        for i in range(self.nr_):
            m_in_r = np.sum(mf[:,i])
            m_out_r = np.sum(mf[i,:])
            m_in_r_tot = m_in_r + m_in[i]
            m_out_r_tot = m_out_r + m_out[i]
            leps[i] = (m_in_r_tot - m_out_r_tot) / max(m_out_r_tot, m_in_r_tot)
            if leps[i] > tol:
                raise ValueError(f"Reactor {i} has an imbalance of {leps[i]}")
            
        # Add attributes
        self.geps_ = geps
        self.leps_ = leps

        # Print information
        if verbose:
            print("Global imbalance of the network: ", self.geps_)
            print("Max imbalance between reactors: ", max(self.leps_))

        return self

    def getReactorAttributes(self, kincorr=False, corrmethod='mean', 
                       verbose=False):
        
        """
        Compute and set the attributes for reactors including volumes, masses, and temperatures.

        Parameters:
        -----------
        kincorr (bool): Flag to enable kinetic corrections.
        corrmethod (str): Method for kinetic corrections ('mean' or 'max').
        verbose (bool): Flag to enable verbose output.

        Returns:
        --------
        self: The object instance with updated reactor attributes.
        """

        # Compute the volumes, mass and average temperatures of the reactors
        vr = []     # volume list (m3)
        mr = []     # mass of reactors (kg)
        Tr = []     # Temperature of reactors (K)

        # If kinetic corrections are needed
        if kincorr:
            Tvar = []
            Tmean = []

        for i in range(self.nr_):
            # the volume is the sum of the cells volume
            vvec = self.datadictionary_['volume']['values'][self.labels_==i]
            vr.append(np.sum(vvec))
            # correction for axial symmetric case
            if self.dim_ == 2:
                vr[i] = vr[i] * 2 * np.pi   

            # compute density as averaged density
            rho = np.mean(self.datadictionary_['density']['values'][self.labels_==i])

            # Get mass of reactor
            mr.append(vr[i]*rho)

            # Compute average temperature
            Tvec = self.datadictionary_['solution']['values'][self.labels_==i,0]
            total_weighted_temp = sum(t * v for t, v in zip(Tvec, vvec))
            total_volume = sum(vvec)
            Tri = total_weighted_temp / total_volume
            Tr.append(Tri[0])

            # Now get quantities for kinetic correction if required
            if kincorr:

                # Get temperature variance in the cluster
                Tvarvec = self.datadictionary_['Tvar']['values'][self.labels_==i]

                # If mean for kinetic corrections is
                if corrmethod == 'mean':
                    total_weighted_temp = sum(t * v for t, v in zip(Tvarvec, vvec))
                    Tvari = total_weighted_temp/total_volume
                    Tvar.append(Tvari[0])
                    Tmean.append(Tr[i])

                # If max for kinetic corrections is used
                elif corrmethod == 'max':
                    Tvar.append(np.max(Tvarvec))
                    Tmean.append(np.max(Tr[i]))

        # Print some info 
        if verbose:
            print(vr)
            print(Tr)
            if kincorr:
                print(Tvar)

        # Get mass flowrates in the reactors
        Mf_in = []
        for i in range(self.nr_):
            Mf_in.append(np.sum(self.Mf_[i,:] + self.inlets_mf_[i]))
        if verbose:
            print("Mass flowrates in the reactors")
            print(Mf_in)

        # Update attributes
        self.Tr_ = Tr
        self.vr_ = vr
        self.mr_ = mr
        self.Mf_in_ = Mf_in
        if kincorr:
            self.kincorr_ = True
            self.Tvar_  = Tvar
            self.Tmean_ = Tmean
        else:
            self.kincorr_ = False 

        return self
    
    def setNetsmokeCRN(self, kinfile, thermofile, canteramech,
                       isothermal=True, verbose=False):
        
        """
        Set up a Chemical Reactor Network (CRN) using PyNetsmoke.

        Parameters:
        -----------
        kinfile (str): Path to the kinetics file.
        thermofile (str): Path to the thermodynamics file.
        canteramech (str): Name of the Cantera mechanism.
        isothermal (bool): Flag indicating whether the system is isothermal.
        verbose (bool): Flag to enable verbose output.

        Returns:
        --------
        crn (ReactorNetwork): The created Chemical Reactor Network object.
        self: The object instance with updated attributes.
        """

        # Create reactors for inlets
        n_in = len(np.where(self.inlets_mf_>0)[0])

        # Now create each reactor and update the list
        Tref = 300
        Pref = self.streams_[0]['gas'].P
        Xref = self.streams_[0]['gas'].X
        
        # Create new connections
        new_connections = np.zeros((n_in, 3))
        count = 0

        # Initialize reactors' list
        R_inlets = []

        for i in range(self.nr_):
            if self.inlets_mf_[i] != 0:

                if verbose:
                    print(f"Reactor {i} is an inlet")

                # Update counter
                count += 1

                # Get associated gas object from the inlet_mixtures
                gas = self.inlet_mixtures_[i]

                # Create the inlet reactor
                rin = Reactor(Rtype='PSR', isothermal=isothermal, tau=1e-6, 
                        Mf=self.inlets_mf_[i], P=Pref, 
                        InletMixture=gas, sp_threshold=1e-5,
                        CanteraMech=canteramech, 
                        isinput=True)
        
                # Append reactor to reactor list
                R_inlets.append(rin)
        
                # Update new connections
                new_connections[count-1, 0] = self.nr_ + count - 1 
                new_connections[count-1, 1] = i
                new_connections[count-1, 2] = self.inlets_mf_[i]

        # Create the new mass flowrate matrix
        mass_netsmoke = np.zeros((self.nr_ + n_in, self.nr_ + n_in))
        for i in range(self.nr_):
            for j in range(self.nr_):
                mass_netsmoke[i, j] = self.Mf_[i, j]

        # Then add the mass flowrates from the fake reactors
        for i in range(n_in):
            mass_netsmoke[int(new_connections[i, 0]), int(new_connections[i, 1])] = new_connections[i, 2]

        # Get new internal mass flowrates
        Mf_in_netsmoke = np.zeros(self.nr_ + n_in)
        for i in range(self.nr_+n_in):
            Mf_in_netsmoke[i] = np.sum(mass_netsmoke[:,i])

        # Initialize reactors list
        R_list = []

        # Scan through the reactors and create the reactor list
        for i in range(self.nr_):
            # Check if is an outlet
            isoutlet = False
            if self.outlets_mf_[i] > 0:
                isoutlet = True
                if verbose:
                    print(f"Reactor {i} is an outlet")
            
            # Get volume, temperature and mass flowrate
            vr = self.vr_[i]
            Tr = self.Tr_[i]
            Mf = Mf_in_netsmoke[i]

            # Security check for kinetic corrections
            if self.kincorr_ == True:
                Tmean = self.Tmean_[i]
                Tvar = self.Tvar_[i]
                if Tmean < 1000:
                    kincorr = False
                else:
                    kincorr = True
            else:
                Tmean = None
                Tvar = None
                kincorr = False

            # Set gas object
            gas = ct.Solution(canteramech)
            gas.TPX = Tr, Pref, Xref

            # Create reactor
            r0 = Reactor('PSR', isothermal=isothermal, 
                volume=vr, Mf=Mf, P=Pref, 
                 InletMixture=gas, InitialStatus=gas, sp_threshold=1e-5,
                 CanteraMech=canteramech, isoutput=isoutlet, KinCorr=kincorr,
                 Tmean=Tmean, Tvar=Tvar)
            
            R_list.append(r0)

        # Add the inlet reactors
        for i in range(n_in):
            R_list.append(R_inlets[i])

        # Create the reactor network object
        crn = ReactorNetwork(R_list,
                mass_netsmoke, kinfile, thermofile)
        
        # Write inputs
        crn.WriteNetworkInput()

        # Update mass flowrates
        self.mass_netsmoke_ = mass_netsmoke

        return crn, self

        













        
            


                    



        






            







    

    




                

            




    


    





        



            

            



            
    

