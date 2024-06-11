import numpy as np
import pandas as pd
import networkx as nx
import time

# Generate graph
def GenGraph(filepath=None, verbose=False, option='fluent',
                      connectivities=None, massflowrates=None):
    
    '''This function will generate the graph of the mesh. 
    Note that this function must be used on the file exported using 
    the Fluent UDF contained in the Fluent_UDF folder. Otherwise, you need to specify 
    manually the connectivities and the exchanged mass flowrates (eventually)
    
    Parameters
    ----------
    datapath : (str)
        path to the connectivity file
        
    Option
    ------
    verbose : (bool)
        display info (default = True)

    option : (str)
        option to handle the data. Available options are 'fluent' or 'custom'
        (default =fluent)

    connectivities : (2d array)
        array containing the connectivities between cell ids, i.e. [2 3; 4 5] and so on.
        It will only be considered if custom is specified as option (default = None)

    massflowrates : (1d array)
        array containing the mass flowrates exchanged by the connectivities, i.e. 3e-4
        means that cells 2 3 are exchanging 3e-4 kg/s (w.r.t. the connectivies array) (default = None)

    Returns
    -------
    G: (networkX graph)
        undirected graph of the mesh (weighted if mass flowrates are specified)

    '''

    # Check for option entry
    if option != 'fluent' and option != 'custom':
        raise ValueError("Options specified is not valid. Specify 'fluent' or 'custom'")
    
    if option == 'fluent':

        if filepath == None:
            raise ValueError("Filepath not specified")

        # Extract connectivity file
        dt = pd.read_csv(filepath, sep='\t')
        id0 = dt['c_id0'].values
        id1 = dt['c_id1'].values
        massflowrates = dt['mass_flowrate'].values[id1 != -1]

        # Remove -1 cells (they are boundary cells and connections)
        id0 = id0[id1 != -1].reshape(-1,1)
        id1 = id1[id1 != -1].reshape(-1,1)

        if len(id0) != len(id1) or len(id0) != len(massflowrates):
            raise ValueError("Connectivity file not read correctly")
        
        # Generate connectivity matrix
        connectivities = np.concatenate((id0, id1), axis=1)

    if verbose:
        print("Mesh has ", len(massflowrates), ' connections')
        print("The maximum exchanged mass flowrate is ", max(massflowrates))

    # Initialize graph
    G = nx.Graph()
    if massflowrates is None:
        for edge in connectivities:
            G.add_edge(edge[0], edge[1])
    else:
        i=0
        for edge in connectivities:
            G.add_edge(edge[0], edge[1], weight=massflowrates[i])
            i+=1

    # Check if the graph is connected
    components = list(nx.connected_components(G))
    if len(components) > 1:
        raise Warning("Careful, your mesh graph is disconnected in ", len(components), ' components!!!')

    return G

# Check cluster connectivities
def CheckClusterConnections(G, labels, datadictionary, verbose=False):

    """
    Check the connectivity of a clustering solution given a graph and labels.

    Parameters
    ----------
    G : networkx.Graph
        The graph representing the mesh.
    labels : array-like
        The labels indicating the cluster to which each node belongs.
    datadictionary : dict
        A dictionary containing data information, including cell IDs.

    Option
    ------
    verbose : (bool)
        display info (default = True)

    Returns
    -------
    flag : bool
        True if the clusters are connected, False if there are disconnected components.

    Raises
    ------
    ValueError
        If the length of nodes in the graph is different from the length of labels.

    Notes
    -----
    This function checks the connectivity of a clustering solution given the graph 
    G and the obtained labels. It returns True if the clusters are connected, 
    and False if there are disconnected components.

    """

    # Check shapes
    if len(G.nodes()) != len(labels):
        raise ValueError("Length of nodes is different from length of labels!")
    # Number of clusters
    k = np.max(labels)+1
    # Init flag
    flag = True
    # Start counter for total unconnected components
    counter = 0
    # Check connectivity by scanning subgraphs of clusters
    for i in range(k):
        # Get nodes ids of the cluster
        nodes = datadictionary['cell_ids'][labels==i]
        # Extract the subgraph
        subG = G.subgraph(nodes).copy()
        # Check connected components
        comps = list(nx.connected_components(subG))
        ncs = len(comps)
        if ncs > 1:
            flag = False
            counter = counter + (ncs-1)
            if verbose:
                print("Cluster ", i, " is disconnected in ", ncs, ' components')

    if verbose:
        print("There are ", counter, " unconnected components in total")

    return flag

def ReassignNodes(G, labels, datadictionary, method='Ncells', itmax = 1000,
                  ncells=1000, v_thresh=1e-6, verbose=False):
    
    """
    Reassign nodes in a graph to ensure cluster connectivity using a specified method.
    
    Parameters
    ----------
    G : networkx.Graph
        Graph representing the connectivity of nodes.
    labels : array-like
        Array of cluster labels for each node.
    datadictionary : dict
        Dictionary containing cell IDs and optionally volumes.
    method : str, optional
        Method for reassignment ('Ncells' or 'volume'). Default is 'Ncells'.
    itmax : int, optional
        Maximum number of iterations. Default is 1000.
    ncells : int, optional
        Threshold for the number of cells in a cluster. Default is 1000.
    v_thresh : float, optional
        Volume threshold for volume-based reassignment. Default is 1e-6.
    verbose : bool, optional
        Flag for verbose output. Default is False.
    
    Returns
    -------
    array-like
        Updated labels ensuring cluster connectivity.
    
    Raises
    ------
    ValueError
        If volume-based reassignment method is specified but 'volume' is not a field in datadictionary.
    
    Notes
    -----
    This function reassigns nodes to ensure the clusters in the graph are connected based on the specified method.
    """

     # Get number of clusters
    k = int(np.max(labels)+1)

    # Check initial connectivity
    flag = CheckClusterConnections(G, labels, datadictionary, verbose=True)

    # Check consistency for specified method
    if method == 'volume' or method == 'Volume':
        if 'volume' not in datadictionary:
            raise ValueError("You speficied volume reassignment method, but volume was not a field in the data!")
        else:
            volumes = datadictionary['volume']['values']
            vtot = np.sum(volumes)
            c_volumes = []
            for i in range(k):
                c_volumes.append(np.sum(volumes[labels==i]))

    # If True, exit the function
    if flag == True:
        labels_new = np.copy(labels)
        print("The clusters are already connected!")
        return labels_new
    
    else:

        # Set counters for routine
        it = 1; counter = 0

        # Initialize convergence flag and new number of clusters
        k_new = k
        labels_new = np.copy(labels)
        convergence = False
        
        # Start reassigning cycle
        while it < itmax and convergence == False:

            # If volume method is specified, compute clusters volume
            if method == 'volume' or method == 'Volume':
                c_volumes = []
                for i in range(k_new):
                    c_volumes.append(np.sum(volumes[labels_new==i]))
                
            # Scan through each cluster
            for i in range(k_new):
                cell_clust_i = datadictionary['cell_ids'][labels_new==i]

                # Check if it's empty
                if cell_clust_i is None:
                    mess = f"Cluster {i} is empty, it will be deleted"; print(mess)
                    # If the cluster is empty, update the idx vector
                    for j in range(len(labels_new)):
                        if labels_new[i] > j:
                            labels_new[i] -= 1

                    # Update number of clusters
                    k_new = k - 1
                    break

                # If the cluster is not empty
                else:
                    # Build the subgraph
                    H = G.subgraph(cell_clust_i).copy()

                    # Get the disconnected components
                    bins = list(nx.connected_components(H))
                    binsize = [len(b) for b in bins]

                    # Sort the disconnected graphs according to their size
                    binsize, ind = zip(*sorted(zip(binsize, range(len(binsize))), reverse=True))
                    bins = [list(b) for b in bins]
                    nbins = len(bins)

                    # If nbins > 1, scan through the subgraph components
                    if nbins > 1:
                        if binsize[1] > 1:
                            # Get cell ids in the second component of the subgraph (the first is the larger)
                            cell_comp = np.array(bins[1]).astype(int) - 1   # Be consistent with Python
                            # Initialize flag for reassignment
                            reassign = False
                            # Reassignment method based on the number of cells
                            if method == 'Ncells':
                                ncells_i = len(cell_comp)
                                if ncells_i < ncells:
                                    reassign = True

                            # Reassignment method based on volume
                            elif method == 'volume' or method == 'Volume':
                                # Compute volume of the component
                                vol_i = np.sum(datadictionary['volume']['values'][cell_comp])
                                # Check if it is larger than volume threshold
                                if vol_i < v_thresh * vtot:
                                    reassign = True

                            # Reassign the cluster if flag is true
                            if reassign:
                                if verbose:
                                    print("Reassigning component of subgraph of cluster ", str(i))

                                # Check for external clusters
                                ext_clust = []
                                for l in range(len(cell_comp)):
                                    neighbour_cells = list(G.neighbors(bins[1][l]))         # Ids of cells (consistent with Fluent)
                                    neighb_ids = np.array(neighbour_cells).astype(int)-1    # Ids of cells to look in labels (consistent with Python)
                                    # Neighbour clusters
                                    neighb_clust = labels_new[neighb_ids]
                                    # Check if any element of neighb_clust is not equal to i (external cluster)
                                    for c in neighb_clust:
                                        if c != i:
                                            ext_clust.append(c)

                                # Check for unique values and counts of groups
                                unique, counts = np.unique(ext_clust, return_counts=True)
                                # Select the cluster that shares more connections
                                ii = np.argmax(counts); best_neighb = unique[ii]
                                # Reassign cells to that cluster
                                labels_new[cell_comp] = best_neighb
                                if verbose:
                                    print("Component was reassigned to cluster ", best_neighb)
                                break
                                
                            # In this case we create a new cluster
                            else:
                                # Update number of clusters
                                k_new = k_new + 1
                                labels_new[cell_comp] = k_new
                                if verbose:
                                    print("A new cluster has been created")
                                break

                        else:
                            if verbose:
                                print("Cluster with only one cell found")
                            cell_i = bins[1][0]
                            neighb = list(G.neighbors(cell_i))
                            neighb = np.array(neighb).astype(int) - 1
                            clust_neighb = labels_new[neighb]
                            new_clust = clust_neighb[0]
                            # Update cluster
                            labels_new[int(cell_i)-1] = new_clust
                            k_new = np.max(labels_new)
                            break 

            # Check the new connectivity
            flag = CheckClusterConnections(G, labels_new, datadictionary, verbose=False)
            if flag == True:
                convergence = True
                print("Convergence reached, the clusters are now connected")
            
            else:
                if it > itmax:
                    convergence = True
                    print("Convergence reached because of iteration limit")
                else:
                    it += 1

    return labels_new

def MatchNClusters(G, labels, datadictionary, nk, method='volume', itmax=1000, verbose=False):

    """
    Adjust the number of clusters in a graph to match a specified number (nk) using either volume-based 
    reassignment or another specified method. The function supports both agglomeration (merging clusters) 
    and splitting clusters to achieve the desired number of clusters.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    labels : array-like
        Array of cluster labels for the nodes in the graph.
    datadictionary : dict
        Dictionary containing additional data, such as 'volume' of nodes.
    nk : int
        Desired number of clusters.
    method : str, optional
        Method for reassignment, default is 'volume'.
    itmax : int, optional
        Maximum number of iterations allowed.

    Returns
    -------
    array-like
        Array of new cluster labels after adjustment.

    Raises
    ------
    ValueError
        If volume-based reassignment is specified but 'volume' is not a field in datadictionary.
    Warning
        If the number of maximum iterations is reached.

    Notes
    -----
    This function adjusts the clustering solution to match the desired number of clusters by either 
    merging or splitting clusters. The volume-based reassignment method requires a 'volume' field in the datadictionary.
    """

    # Preliminary checks
    if method=='volume' or method=='Volume':
        if 'volume' not in datadictionary:
            raise ValueError("Volume-based reassignment was specified, but volume is not a field in datadictionary!")

    # Get current number of clusters
    k = np.max(labels)+1

    # Check situation
    if k == nk:
        print("Number of clusters matches number of reactors, exiting the function")
        return labels
    elif k > nk:
        print("Number of clusters greater than desired number of clusters, agglomeration is needed")
    else:
        print("Number of clusters lower than desired number of clusters, splitting is needed")

    # Remove empty clusters
    labels_new = RemoveEmptyClusters(labels)
    k_new = np.max(labels_new)+1

    # If agglomeration is needed
    if k > nk:
        # Start iteration counter 
        iter = 0
        if verbose:
            print("--- Reactor matching code iteration ", iter, ' ---')
        # While cycle until k matches nk
        while k_new > nk and iter < itmax:
            
            # Calculate volumes of the clusters
            cvol = []
            for i in range(k_new):
                vcells_i = datadictionary['volume']['values'][labels_new==i]
                cvol.append(np.sum(vcells_i))
            
            # Find smallest cluster
            sid = np.argmin(cvol)
            vmin = np.min(cvol)
            if vmin == 0:
                print(f"Cluster {sid} has volume 0, must be removed")


            # Get the ids of the cells contained in sid cluster
            cell_ids = datadictionary['cell_ids'][labels_new==sid]       # Ids of cells consistent with Fluent
            
            # Check for external clusters
            ext_clust = []
            for l in range(len(cell_ids)):
                # Get neighbor cells
                neighbour_cells = list(G.neighbors(cell_ids[l]))            # Ids of cells (consistent with Fluent)
                neighb_ids = np.array(neighbour_cells).astype(int)-1        # Ids of cells to look in labels (consistent with Python)
                # Neighbour clusters
                neighb_clust = labels_new[neighb_ids]
                # Check if any element of neighb_clust is not equal to i (external cluster)
                for c in neighb_clust:
                    if c != sid:
                        ext_clust.append(c)

            # Check for unique values and counts of groups
            unique, counts = np.unique(ext_clust, return_counts=True)
            if len(unique) == 0:
                raise ValueError(f"Cluster {i} does not have any neighbour. Not possible?")
            # Select the cluster that shares more connections
            ii = np.argmax(counts); best_neighb = unique[ii]

            # Reassign cells to that cluster
            labels_new[labels_new==sid] = best_neighb

            # Update all indices greater than sid 
            labels_new[labels_new > sid] = labels_new[labels_new > sid] - 1

            # Update k_new
            k_new = np.max(labels_new) + 1

            # Update iteration counter
            iter = iter + 1
            if iter == itmax:
                raise Warning("Number of maximum iterations reached")
            
    # Check connection
    flag = CheckClusterConnections(G, labels_new, datadictionary, verbose=verbose)
    if flag == False:
        raise ValueError("Now the clustering solution is not connected")
            
    print("Actual number of clusters = ", k_new)
    return labels_new
            
def RemoveEmptyClusters(labels):
    # Determine the range of possible cluster values
    max_cluster = np.max(labels)
    all_clusters = set(range(max_cluster + 1))
    
    # Find the unique clusters present in the labels
    present_clusters = set(np.unique(labels))
    
    # Identify the missing clusters
    missing_clusters = sorted(all_clusters - present_clusters)
    
    # Print the missing clusters
    if missing_clusters:
        print("Removed clusters:", missing_clusters)
    else:
        print("No clusters removed")

    # Create a mapping from old cluster numbers to new cluster numbers
    cluster_mapping = {old_cluster: new_cluster for new_cluster, old_cluster in enumerate(present_clusters)}
    
    # Apply the mapping to the labels array
    new_labels = np.array([cluster_mapping[label] for label in labels])
    
    return new_labels




    











                                

                            
                                
                                                                






                        




                    




                        






                




    
