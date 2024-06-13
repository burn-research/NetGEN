import numpy as np
import pandas as pd
import cantera as ct

def GenDataFluent(datainfo, datapath, dim, verbose=False):

    '''This function will import the data from the selected CFD files contained 
    in the dictionary datainfo. The datafiles are .csv or ascii files stored
    in the datapath folder. Note the this function must be used only to handle
    files exported from Fluent

    Parameters
    ----------
    datainfo : (dictionary)
        example datainfo['fieldname'] = 'volume'
                datainfo['filename'] = 'volumefile.csv
    datapath : (str)
        path to where the data are stored
    dim : (int)
        dimension of the geometrical CFD (2 for 2D and 3 for 3D)

    Returns
    -------
    datadictionary : (dictionary)
        dictionary with also numerical values'''
    
    # Check dim entry
    if dim != 2 and dim != 3:
        raise ValueError("Dimension incorrect. Specify 2 for 2D or 3 for 3D")
    # Index to start exporting numerical values from columns
    istart = dim+1
    # Initialize the dictionary of the data
    datadictionary = {}
    # Export all the data contained in the files
    i=0
    for key in datainfo:
        filename = datapath + "/" + datainfo[key]
        dt = pd.read_csv(filename, sep=',')
        datadictionary[key] = {}
        # Create two fields for each key
        datadictionary[key]['values'] = dt.values[:,istart:]

        # Adjust names by removing initial space left by Fluent
        names = list(dt.columns[istart:])
        newnames = []
        for name in names:
            newnames.append(name.split()[0])

        datadictionary[key]['names']  = newnames
        # Create one voice for coordinates
        if i == 0:
            datadictionary['coordinates'] = {}
            datadictionary['coordinates']['values'] = dt.values[:,1:istart]
            datadictionary['coordinates']['names'] = dt.columns[istart:]
            datadictionary['cell_ids'] = dt.values[:,0]
            i+=1

    if verbose:
        print("Keys present in the datadictionary are:")
        for key in datadictionary:
            print(key, ", shape = ", np.shape(datadictionary[key]['values']))
        print("")
        print("Names of solution data are:")
        print(datadictionary['solution']['names'])
    
    return datadictionary

# ------------- Utilities functions to calculate useful quantities ------------ #
def getMixtureFraction(datadictionary, fuel, oxidizer, canteramech="gri30.yaml", basis="mole", datainterface="fluent"):

    # Check if datadictionary contains the solution field
    if "solution" not in datadictionary:
        raise ValueError("The 'solution' field is not in datadictionary. Please, provide a thermo-chemical solution before extracting mixture fraction")

    # Initialize cantera object
    gas = ct.Solution(canteramech)

    # Get composition
    comp = datadictionary['solution']['values'][:,1:]

    # Initialize mixture fraction
    npts = np.shape(comp)[0]
    Z = np.zeros(npts)

    # Get column names according to cantera with fluent data interface
    if datainterface == "fluent":
        # Names of columns 
        names = datadictionary['solution']['names'][1:]
        # Initialize list of Fluent names
        fnames = []
        # Convert to capital and species names
        for name in names:
            if '-' in name:
                fnames.append(name.split('-')[1].upper())
            else:
                fnames.append(name.upper())
        
        # Iterate and calculate point-wise mixture fraction
        for i in range(npts):
            x = np.zeros(gas.n_species)
            for j in range(len(fnames)):
                idsp = gas.species_index(fnames[j])
                x[idsp] = comp[i,j]
            if basis == "mole":
                gas.TPX = 300.0, 101325.0, x
                # Calculate mixture fraction
                z = gas.mixture_fraction(fuel, oxidizer, basis=basis, element="Bilger")
            else:
                gas.TPY = 300.0, 101325.0, x
                # Calculate mixture fraction
                z = gas.mixture_fraction(fuel, oxidizer, basis=basis, element="Bilger")
            # Update mixture fraction array
            Z[i] = z

    else:
        raise ValueError("datainterface not recognized")
    
    return Z



        


            




        

    



        




    


    


    




    
