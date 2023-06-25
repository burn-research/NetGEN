# NetGEN
*NetGEN* is a Matlab suite for automatically generating Chemical Reactor Networks (CRNs) from CFD data. 
The package is based on connectivity-enhanced unsupervised learning clustering algorithms. In particular, 
clustering algorithms are used to post-process the data stored in the computational cells of a CFD simulation,
finding zones with similar thermo-chemical features that can be modelled as a series of interconnected chemical reactors.
To guarantee that the identified zones represent geometrically connected areas within the combustor, graph algorithms
are employed to reassign the cells, making sure that each cluster is composed of spatially connected zones. 
This will guarantee the conservation of mass in the CRN model. For more details, see the publication [[1]](#1).

<a id="1">[1]</a> 
M. Savarese, A. Cuoci, W. De Paepe, A. Parente (2023)
Machine Learning clustering algorithms for the automatic generation of chemical reactor networks from CFD simulations
Fuel, 343, 127945


## Main features of the code
Several functionalities are available in the provided code, in particular:
- Functions for the pre-processing of CFD data are provided
- A FLUENT interface is available for extracting the required files via UDF functions
- Unsupervised Learning algorithms and data pre-processing tools are available (k-means and VQPCA)
- Graph algorithms to perform connectivity scanning and nodes reassignment are provided and can be of more general purpose
- An interface with the CRN solver **NetSMOKE++** is provided. The solver can be available upon request
- Post-processing tools to plot and visualize CRN results

A detailed description of the tools is available in the **Example** folder, where a **Matlab** live script is provided

## Usage and License
The code is open source. Please cite [[1]](#1).


