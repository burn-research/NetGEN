# NetGEN

**Automatic Generation of Chemical Reactor Networks from CFD data**

**NetGEN** is a suite for automatically generating Chemical Reactor Networks (CRNs) from CFD data. 
The package is based on connectivity-enhanced unsupervised learning clustering algorithms. In particular, 
clustering algorithms are used to post-process the data stored in the computational cells of a CFD simulation,
finding zones with similar thermo-chemical features that can be modelled as a series of interconnected chemical reactors.
To guarantee that the identified zones represent geometrically connected areas within the combustor, graph algorithms
are employed to reassign the cells, making sure that each cluster is composed of spatially connected zones. 
This will guarantee the conservation of mass in the CRN model. For more details, see the publication [[1]](#1).

## Table of Contents

- [Python Installation](#pythoninstallation)
- [Matlab Installation](#matlabinstallation)
- [Usage](#usage)
- [Features](#features)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Python Installation
1. Clone the repository:
    ```sh
    git clone https://github.com/burn-research/NetGEN.git
    ```

2. Enter the Python source code folder:
    ```sh
    cd Python
    ```

3. Install:
    ```sh
    python setup.py install
    ```

## Matlab Installation
Open Matlab, then type:
```sh
addpath(genpath('MATLAB/src/'))
```

## Usage
The code is directly usable. To run examples, you can either find examples in the
Python/examples/ folder or in the MATLAB/Example/ folder. Those example contains
a series of options that can be customized.

## Features
- Unsupervised Learning (k-means, vector quantization PCA)
- Graph-based clustering
- Fluent - NetSMOKE++ interface for CRN

## Acknowledgements
We want to deeply thanks the Creck group at Politecnico di Milano for their collaboration 
to this project, in particular for the development of the CRN solver NetSMOKE++. 
The CRN solver is available upon request for academic usage

## License
The code is open source. Please cite [[1]](#1).

<a id="1">[1]</a> 
M. Savarese, A. Cuoci, W. De Paepe, A. Parente (2023)
Machine Learning clustering algorithms for the automatic generation of chemical reactor networks from CFD simulations
Fuel, 343, 127945


