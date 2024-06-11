import numpy as np
import pandas as pd

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

        

    



        




    


    


    




    
