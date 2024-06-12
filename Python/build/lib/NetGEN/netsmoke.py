import numpy as np
import cantera as ct
import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import time

class Reactor:
    def __init__(self, Rtype, isothermal=False, volume=None, tau=None, Mf=None, P=None, 
                 L=None, D=None, Ac=None, Aq=None, Uq=None, Tenv=None, AcOverP=None,
                 InletMixture=None, InitialStatus=None, sp_threshold=1e-5,
                 isinput=False, isoutput=False,
                 CanteraMech='gri30.cti', KinPath="dummy",
                 KinFile=None, ThermoFile=None, PreProcessor=False,
                 KinCorr=False, Tmean=None, Tvar=None, CorrType="beta"):
        
        """
        Initialize a Reactor object with specified parameters. Valid for OpenSMOKEpp and NetSMOKEpp

        Parameters:
        -----------
        Rtype (str): Type of the reactor ('PSR' or 'PFR').
        isothermal (bool): Flag indicating whether the reactor is isothermal.
        volume (float): Reactor volume [m3].
        tau (float): Reactor residence time [s].
        Mf (float): Mass flowrate [kg/s].
        P (float): Pressure [Pa].
        L (float): Length [m].
        D (float): Diameter [m].
        Ac (float): Cross-sectional area [m2].
        Aq (float): Heat exchange area [m2].
        Uq (float): Heat transfer coefficient [W/m2/K].
        Tenv (float): Environment temperature [K].
        AcOverP (float): Cross section over perimeter [m].
        InletMixture: Inlet mixture (Cantera object).
        InitialStatus: Initial status (Cantera object).
        sp_threshold (float): Threshold to write initial composition dictionary.
        isinput (bool): True if the reactor is an inlet reactor.
        isoutput (bool): True if the reactor is an outlet reactor.
        CanteraMech (str): Chemical mechanism file for Cantera.
        KinPath (str): Kinetics folder (for pre-processed kinetic mechanism).
        KinFile (str): Path to kinetic mechanism file (usually .dat, .kin...).
        ThermoFile (str): Path to thermodynamic file.
        PreProcessor (bool): Flag indicating if kinetic and thermo file are specified.
        KinCorr (bool): Flag for kinetic corrections activation.
        Tmean (float): Average temperature for kinetic corrections.
        Tvar (float): Average variance for kinetic corrections.
        CorrType (str): Type of kinetic correction to apply (default 'beta').

        Raises:
        -------
        ValueError: If the reactor type is not 'PSR' or 'PFR'.
        """
        
        # Check if Rtype exists
        if Rtype != 'PSR' and Rtype != 'PFR':
            raise ValueError('Unknown reactor type. Specify PSR or PFR (must be specified upper case)')
        
        self.Rtype = Rtype                          # Reactor type ('PSR or PFR)
        self.isothermal = isothermal                # Isothermal or Non-isothermal 
        # Reactor parameters
        self.volume = volume                        # Reactor volume            [m3]
        self.tau = tau                              # Reactor residence time    [s]
        self.Mf = Mf                                # Mass flowrate             [kg/s]
        self.P = P                                  # Pressure                  [Pa]
        self.L = L                                  # Length                    [m]
        self.D = D                                  # Diameter                  [m]
        self.Ac = Ac                                # Cross sectional area      [m2]
        self.Aq = Aq                                # Heat exchange area        [m2]
        self.Uq = Uq                                # Heat transfer coefficient [W/m2/K]
        self.Tenv = Tenv                            # Environment temperature   [K]
        self.AcOverP = AcOverP                      # Cross section over per.   [m]
        # Reactor initialization and inlets
        self.InletMixture = InletMixture            # Inlet mixture (cantera object)
        self.InitialStatus = InitialStatus          # InitialStatus (cantera object)
        self.sp_threshold = sp_threshold            # Threshold to write initial composition dictionary
        self.CanteraMech = CanteraMech              # Chemical mechanism file cantera
        # Only for reactor networks
        self.isinput = isinput                      # True if is an inlet reactor
        self.isoutput = isoutput                    # True if is an outlet reactor
        # ROPA options
        self.ropa = False                           # True when you do SetROPA
        self.RopaThreshold = 0.001                  # Relative threshold for reaction rate
        self.RopaSpecies = None                     # ROPA species
        self.RopaThreshold = None                   # Threshold for ROPA species
        # Kinetic mechanism informations
        self.KinPath = KinPath                      # Kinetics Folder (to pre-processed kinetic mechanism)
        self.KinFile = KinFile                      # Path to kinetic mechanism file (usually .dat, .kin...)
        self.ThermoFile = ThermoFile                # Path to thermodynamic file
        self.PreProcessor = PreProcessor            # If kinetic and thermo file are specified, it must be true
        # Kinetic corrections
        self.KinCorr  = KinCorr                     # Bool for kinetic corrections actvation
        self.Tmean    = Tmean                       # Average temperature for kinetic corrections
        self.Tvar     = Tvar                        # Average variance for kinetic corrections
        self.CorrType = CorrType                    # Type of kinetic correction to apply (default beta)

        # Initialize file path
        self.FilePath = None


    def SetROPA(self, species, reference, threshold=0.001):
        """
        Set up a Rate of Production Analysis (ROPA).

        Parameters:
        -----------
        species (str): Species for ROPA analysis.
        reference (str): Reference species for ROPA analysis.
        threshold (float): Relative threshold for reaction rate.

        Returns:
        --------
        self: The Reactor object with ROPA settings.
        """

        # Activate flag
        self.ropa = True
        self.RopaSpecies = species
        self.ReferenceSpecies = reference
        self.RopaThreshold = threshold
        return self

    def WriteInput(self, filepath=None, Rid=None):

        '''
        Write the reactor input dictionary file.

        Parameters:
        -----------
        filepath (str, optional): Path where the input file will be saved.
        Rid (int, optional): Reactor ID used for naming the input file.

        Returns:
        --------
        self: The Reactor object itself.

        Raises:
        -------
        ValueError: If the path to thermo or kinetic mechanism file is not specified
                    when the preprocessor is selected.
        '''

        # Check if Rid was specified
        if Rid != None:
            self.Rid = Rid
        else:
            self.Rid = None
        # Check if filepath was specified, otherwise save it where you are
        if filepath != None:
            self.FilePath = filepath
        else:
            self.FilePath = '.'
            
        # PSR input dictionary
        if self.Rtype == 'PSR':
            if filepath == None and Rid == None:
                filename = 'input.dic'
            # If Rid was specified, change name to input file
            elif filepath == None and isinstance(Rid, int):
                filename = 'input.cstr.' + str(Rid) + '.dic'
            # If also the path was specified save it there
            elif isinstance(filepath, str) and isinstance(Rid, int):
                filename = filepath + '/input.cstr.' + str(Rid) + '.dic'
            # If path was specified but not Rid, save it there just as input.dic
            elif isinstance(filepath, str) and Rid == None:
                if os.path.exists(filepath) == False:
                    os.mkdir(filepath)
                filename = filepath + '/input.dic'

        # PFR input dictionary (same rules as before)
        elif self.Rtype == 'PFR':
            if filepath == None and Rid == None:
                filename = 'input.dic'
            elif filepath == None and isinstance(Rid, int):
                filename = 'input.pfr.' + str(Rid) + '.dic'
            elif isinstance(filepath, str) and isinstance(Rid, int):
                filename = filepath + '/input.pfr.' + str(Rid) + '.dic'
            elif isinstance(filepath, str) and Rid == None:
                if os.path.exists(filepath) == False:
                    os.mkdir(filepath)
                filename = filepath + '/input.dic'

        # Open file
        with open(filename, 'w') as f:
            # Check if it is a PSR or a PFR
            if self.Rtype == 'PSR':
                f.write('Dictionary PerfectlyStirredReactor \n { \n')
            elif self.Rtype == 'PFR':
                f.write('Dictionary PlugFlowReactor \n { \n')

            # Write attributes
                
            # Options for PSR
            if self.isothermal == True and self.Rtype == 'PSR':
                f.write('     @Type     Isothermal-ConstantPressure;\n')
            elif self.isothermal == False and self.Rtype == 'PSR':
                f.write('     @Type     NonIsothermal-ConstantPressure;\n')
            # Options for PFR
            elif self.isothermal == True and self.Rtype == 'PFR':
                f.write('     @Type     Isothermal;\n')
                f.write('     @ConstantPressure   true;\n')
            elif self.isothermal == False and self.Rtype == 'PFR':
                f.write('     @Type     NonIsothermal;\n')
                f.write('     @ConstantPressure   true;\n')
            else:
                 f.write('     @Type     Isothermal-ConstantPressure;\n')
                # raise ValueError("Unable to specify @Type correctly. Specify if isothermal as bool and the type as string ('PSR' or 'PFR')")

            # Kinetic dictionary
            if self.PreProcessor == False:
                f.write('     @KineticsFolder     %s;\n' %self.KinPath)
            elif self.PreProcessor == True:
                f.write('      @KineticsPreProcessor     kinetic-mechanism; \n')
                if self.KinFile == None or self.ThermoFile == None:
                    raise ValueError('Preprocessor is selected but path to thermo or kinetic mech is not specified!')

            # Value instances
            # Volume
            if isinstance(self.volume, (float, int)):
                f.write('     @Volume    %e m3;\n' %self.volume)
            # Residence time
            if isinstance(self.tau,  (float, int)):
                f.write('     @ResidenceTime    %e s;\n' %self.tau)
            # Mass flowrate
            if isinstance(self.Mf,  (float, int)):
                f.write('     @MassFlowRate    %e kg/s;\n' %self.Mf)
            # Length
            if isinstance(self.L,  (float, int)):
                f.write('     @Length    %e m;\n' %self.L)           
            # Diameter
            if isinstance(self.D,  (float, int)):
                f.write('     @Diameter    %e m;\n' %self.D)
            # Heat exchange area
            if isinstance(self.Aq,  (float, int)):
                f.write('     @ExchangeArea    %e m2;\n' %self.Aq)        
            # Heat transfer coefficient
            if isinstance(self.Uq,  (float, int)):
                f.write('     @GlobalThermalExchangeCoefficient    %e W/m2/K;\n' %self.Uq) 
            # Environment temperature
            if isinstance(self.Tenv,  (float, int)):
                f.write('     @EnvironmentTemperature    %d K;\n' %self.Tenv)
            # Cross section over perimeter
            if isinstance(self.AcOverP,  (float, int)):
                f.write('     @CrossSectionOverPerimeter    %d m;\n' %self.AcOverP)

            # Initial status
            if self.InitialStatus != None and self.isothermal == False:
                f.write('     @InitialStatus     initial-status;\n')
            
            # Inlet status
            if self.InletMixture != None:
                f.write('     @InletStatus     inlet-mixture;\n')

            # ROPA option
            if self.ropa == True:
                f.write('     @OnTheFlyROPA      ropa;\n')

            # Kinetic corrections
            if self.KinCorr == True:
                f.write('     @KineticCorrections      kinetic-corrections;\n')

            # Close dictionary
            f.write('} \n')

            ### Additional dictionaries ####
            if self.InletMixture != None:
                M = self.InletMixture   # Cantera object
                # Extract main quantities
                T = M.T
                P = M.P
                # Extract composition
                X = M.X
                # Extract names
                names = M.species_names
                # Select only a small portion of species to be written
                mfsp = []
                spls = []
                thresh = 1e-6
                for i in range(len(names)):
                    if X[i] > thresh:
                        mfsp.append(X[i])
                        spls.append(names[i])
                f.write('Dictionary inlet-mixture \n { \n')
                f.write('     @Temperature     %e K;\n' % T)
                f.write('     @Pressure        %e Pa;\n' % P)
                f.write('@Moles     ')
                # Write species
                for i in range(len(mfsp)):
                    if i < len(mfsp)-1:
                        f.write('%s %f \n' % (spls[i], mfsp[i]))
                    else:
                        f.write('%s %f; \n' % (spls[i], mfsp[i]))
                f.write('} \n')

            ### Additional dictionaries ####
            if self.InitialStatus != None and self.isothermal == False:
                M = self.InletMixture   # Cantera object
                # Extract main quantities
                T = M.T
                P = M.P
                # Extract composition
                X = M.X
                # Extract names
                names = M.species_names
                # Select only a small portion of species to be written
                mfsp = []
                spls = []
                thresh = 1e-6
                for i in range(len(names)):
                    if X[i] > thresh:
                        mfsp.append(X[i])
                        spls.append(names[i])
                f.write('Dictionary initial-status \n { \n')
                f.write('     @Temperature     %e K;\n' % 2500.0)
                f.write('     @Pressure        %e Pa;\n' % P)
                f.write('     @Moles     ')
                # Write species
                for i in range(len(mfsp)):
                    if i < len(mfsp)-1:
                        f.write('%s %f \n' % (spls[i], mfsp[i]))
                    else:
                        f.write('%s %f; \n' % (spls[i], mfsp[i]))
                f.write('} \n')

            # ROPA dictionary
            if self.ropa == True:
                f.write('Dictionary ropa \n{ \n')
                f.write('     @Species ')
                for i in range(len(self.RopaSpecies)):
                    f.write(' %s' % self.RopaSpecies[i])
                f.write(';\n')
                # Reference species
                f.write('     @ReferenceSpecies ')
                for i in range(len(self.ReferenceSpecies)):
                    f.write(' %s' % self.ReferenceSpecies[i])
                f.write(';\n')
                # Threshold
                f.write('     @Threshold    %e;\n' % self.RopaThreshold)
                f.write('} \n')

            # Kinetic corrections dictionary
            if self.KinCorr == True:
                f.write('\nDictionary kinetic-corrections \n{ \n')
                f.write('     @AverageTemperature     %e;\n' % self.Tmean)
                f.write('     @TemperatureVariance     %e;\n' % self.Tvar)
                f.write('     @Type     %s; \n' % self.CorrType)
                f.write('} \n')

            # Kinetic mechanism dictionary
            if self.PreProcessor == True:
                f.write('\nDictionary kinetic-mechanism \n{ \n')
                f.write('     @Kinetics     %s; \n' % self.KinFile)
                f.write('     @Thermodynamics    %s;   \n' % self.ThermoFile)
                f.write('     @Output     kinetics; \n } \n')

        return self
    
    def RunSimulation(self, ospath):

        '''
        Runs the simulation of the single reactor using the OpenSMOKE++ executable files.

        Parameters:
            ospath (str): The path to the OpenSMOKE++ installation directory.

        Raises:
            ValueError: If the specified path does not exist.

        Returns:
            Reactor: The Reactor object.
        '''

        if os.path.exists(ospath) == False:
            raise ValueError("Your specified path %s does not exists!" % ospath)

        # Generate the command
        if self.Rtype == 'PSR':
            command = 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:' + ospath +'/lib ; ' + ospath + '/bin/OpenSMOKEpp_PerfectlyStirredReactor.sh --input input.dic'
        elif  self.Rtype == 'PFR':
            command = 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:' + ospath +'/lib ; ' + ospath + ' /bin/OpenSMOKEpp_PlugFlowReactor.sh --input input.dic'

        # Logfile
        log_file = "logfile.txt"

        print("Executing NetSMOKEpp...")
        # Execute command
        if self.FilePath != None:
            os.chdir(self.FilePath)
            os.system(f"{command} >> {log_file}")
            os.chdir('../')
        else:
            os.system(f"{command} >> {log_file}")

        return self

    def ExtractRopa(self, filepath=None):

        '''
        Extracts the ROPA (Rate of Production Analysis) output as an array of net reaction rates.

        Parameters:
            filepath (str, optional): The path to the ROPA output file. If not provided and `Rid` is not None, the default path is set to 'ReactorNetwork/Output/Reactor.Rid/ROPA.out'. If `Rid` is None or no `filepath` is provided, it defaults to 'FilePath/Output/ROPA.out'.

        Returns:
            list: An array containing the net reaction rates extracted from the ROPA output file.
        '''

        if filepath == None and self.Rid != None:
            filepath = 'ReactorNetwork/Output/Reactor.' + str(self.Rid) + '/ROPA.out'
        else:
            filepath = self.FilePath + '/Output/ROPA.out'

        net_rates = []
        for sp in self.RopaSpecies:
            with open(filepath, 'r') as f:
                spcounter = 0
                for line in f:
                    ss = line.split()
                    if len(ss) > 0:
                        if ss[0] == sp:
                            spcounter += 1
                            if spcounter == 2:
                                net_rates.append(float(ss[6]))

        return net_rates
    
    def ExtractInletComp(self, filepath=None, sp_threshold=1e-10):

        '''
        Extracts the inlet status of the mixture to a certain reactor in a reactor network.

        This function is specifically designed for reactors in reactor networks.

        Parameters:
            filepath (str, optional): The path to the log file containing inlet status information. 
                If not provided, the default path is set to 'ReactorNetwork/Output/Reactor.Rid/log.inlet', where `Rid` is the reactor ID.

            sp_threshold (float, optional): The threshold value below which species mole 
                fractions are considered negligible. Defaults to 1e-10.

        Returns:
            gas: A Cantera Solution object representing the mixture's composition at the inlet of the specified reactor.
        '''

        if self.Rid == None:
            raise ValueError("This function is only for reactors in reactor networks!")

        if filepath == None:
            filepath = 'ReactorNetwork/Output/Reactor.' + str(self.Rid) + '/log.inlet'
        
        # Extract file
        df = pd.read_csv(filepath, sep='\s+')

        # Export T, P
        T = df.values[-1,1]
        P = 101325.0

        # In columns, locate all the strings that have '_x" inside
        sp_list = []
        Y_val   = []
        for i, col in enumerate(df.columns[2:], 2):
            ss = col.split('(')
            if len(ss) > 1:
                if df.values[-1,i] < 1:
                    sp_list.append(ss[0])
                    Y_val.append(df.values[-1,i])

        # Create solution
        gas = ct.Solution(self.CanteraMech)
        # Extract number of species
        ns = gas.n_species
        # sp_list = gas.species_names
        Y_old = gas.Y
        Y_new = Y_old
        for i, sp in enumerate(sp_list):
            if Y_val[i] > sp_threshold:
                Y_new[gas.species_index(sp)] = Y_val[i]
            else:
                Y_new[i] = 0.0

        # Set new gas state
        gas.TPY = T, P, Y_new

        return gas
    
    def ExtractOutput(self, filepath=None):

        '''
        Extracts the output of a certain reactor and returns it as a Cantera solution object.

        Parameters:
            filepath (str, optional): The path to the output file. 
                If not provided, the default path is set to 'ReactorNetwork/Output/Reactor.Rid/Output.out', 
                where `Rid` is the reactor ID.

        Returns:
            gas: A Cantera Solution object representing the composition and state of the reactor output.
        '''

        if filepath == None and self.Rid != None:
            filepath = 'ReactorNetwork/Output/Reactor.' + str(self.Rid) + '/Output.out'
        else:
            filepath = self.FilePath + '/Output/Output.out'

        df = pd.read_csv(filepath, sep='\s+')

        # Export T, P
        T = df.values[-1,4]
        P = df.values[-1,5]

        # In columns, locate all the strings that have '_x" inside
        sp_list = []
        X_val   = []
        for i, col in enumerate(df.columns):
            ss = col.split('(')
            if len(ss) > 1 and ss[0] != "rho[kg/m3]":
                if df.values[-1,i] < 1:
                    sp_list.append(ss[0])
                    X_val.append(df.values[-1,i])

        # Create solution
        gas = ct.Solution(self.CanteraMech)
        # Extract number of species
        ns = gas.n_species
        sp_list = gas.species_names
        X_old = gas.X
        X_new = X_old
        for i, sp in enumerate(sp_list):
            X_new[gas.species_index(sp)] = X_val[i]

        # Set new gas state
        gas.TPX = T, P, X_new

        return gas
    
### Create class for reactor networks
class ReactorNetwork:

    """
    Class for representing Reactor Networks in NetSMOKE++ 

    Attributes:
        Rlist (list): List of reactors.
        MassFlowrates (ndarray): Matrix of internal mass flowrates.
        KinFile (str): Path to the kinetic file.
        ThermoFile (str): Path to the thermodynamic file.
        MinIt (int): Minimum number of iterations.
        MaxIt (int): Maximum number of iterations.
        AtomicThreshold (float): Atomic error threshold.
        NonIsothermalThreshold (float): Non-isothermal error threshold.
        MaxUnbalance (float): Maximum unbalance internal mass flowrates.
        SpeciesMonitor (list): List of species to monitor.
        CanteraMech (str): Cantera chemical mechanism.
        Nr (int): Number of reactors.

    Raises:
        ValueError: If MassFlowrates matrix is non-square or the number of reactors is not consistent between reactor list and mass flowrates.

    """
    def __init__(self, Rlist, MassFlowrates, KinFile, ThermoFile,
                 MinIt=5, MaxIt=500, AtomicThreshold=1e-3, NonIsothermalThreshold=1e-4,
                 MaxUnbalance=0.05, SpeciesMonitor=['NO','CO'], CanteraMech='gri30.cti'):
        
        # Default initializer
        self.Rlist = Rlist                                          # List of reactors
        self.MassFlowrates = MassFlowrates                          # Matrix of internal mass flowrates
        self.KinFile = KinFile                                      # Path to kinetic file
        self.ThermoFile = ThermoFile                                # Path to thermo file
        self.MinIt = MinIt                                          # Minimum number of iterations
        self.MaxIt = MaxIt                                          # Maximum number of iterations
        self.AtomicThreshold = AtomicThreshold                      # Atomic error threshold
        self.NonIsothermalThreshold = NonIsothermalThreshold        # Non isothermal error threshold
        self.MaxUnbalance = MaxUnbalance                            # Max unbalance internal mass flowrates
        self.SpeciesMonitor = SpeciesMonitor                        # List of species to monitor
        self.CanteraMech = CanteraMech                              # Cantera chemical mechanism

        # Check that Rlist and MassFlowrates have same dimensions
        nr = len(Rlist)
        (nr1,nr2) = np.shape(MassFlowrates)

        if nr1 != nr2:
            raise ValueError('Mass flowrates matrix is non square. Check!')

        if nr != nr1:
            raise ValueError('Number of reactors not consistent between reactor list and mass flowrates')
        
        # Number of reactors
        self.Nr = nr1
        
    # Write input file function
    def WriteNetworkInput(self, FolderName='ReactorNetwork', DicName='input.dic'):

        """
        Write input file function.

        This function will write only the global input.dic file for the 
        reactor network simulation. The other reactors, you must write their input file
        singluarly.

        Args:
            FolderName (str, optional): Folder name. Defaults to 'ReactorNetwork'.
            DicName (str, optional): Dictionary name. Defaults to 'input.dic'.

        Returns:
            self: Instance of the ReactorNetwork class.

        Raises:
            ValueError: If no inlets or outlets detected.
        """
        # Get current working directory
        cwd = os.getcwd()
        self.WorkingDir = cwd

        OutputFile = FolderName + '/' + DicName
        if os.path.exists(FolderName) == False:
            os.mkdir(FolderName)

        # Write single reactors inputs in FolderName
        idr = 0
        for r in self.Rlist:
            r.WriteInput(filepath=FolderName, Rid=idr)
            idr += 1

        # Identify PSRs and PFRs
        id_psrs = []; id_pfrs = []
        for i in range(len(self.Rlist)):
            if self.Rlist[i].Rtype == 'PSR':
                id_psrs.append(i)
            elif self.Rlist[i].Rtype == 'PFR':
                id_pfrs.append(i)

        # Number of reactors
        Nr = np.shape(self.MassFlowrates)[0]
        # Non zero connections
        Nconns = np.count_nonzero(self.MassFlowrates)

        # Count inlets
        InletList = []
        OutletList = []
        for i in range(Nr):
            if self.Rlist[i].isinput == True:
                InletList.append(i)
            elif self.Rlist[i].isoutput == True:
                OutletList.append(i)
        
        if len(InletList) == 0 or len(OutletList) == 0:
            raise ValueError('No inlets or outlets detected. Specify inlet/outlet reactors')

        # Write output file 
        with open(OutputFile, 'w') as f:
            f.write('Dictionary ReactorNetwork \n { \n')
            f.write('@KineticsPreProcessor     kinetic-mechanism;\n')
            f.write('@MinIterations                     %d;\n' % self.MinIt)
            f.write('@MaxIterations                     %d;\n' % self.MaxIt)
            f.write('@AtomicErrorThreshold              %e;\n' % self.AtomicThreshold)
            f.write('@NonIsothermalErrorThreshold       %e;\n' % self.NonIsothermalThreshold)
            f.write('@MaxUnbalance                      %e;\n' % self.MaxUnbalance)
            # Write perfectly stirred reactors
            if len(id_psrs) > 0:
                f.write('@PerfectlyStirredReactors \n')
                for i in range(len(id_psrs)):
                    if i < len(id_psrs)-1:
                        f.write('%d       input.cstr.%d.dic\n' % (id_psrs[i], id_psrs[i]))
                    else:
                        f.write('%d       input.cstr.%d.dic;\n' % (id_psrs[i], id_psrs[i]))
            # Write plug flow reactors
            if len(id_pfrs) > 0:
                f.write('@PlugFlowReactors \n')
                for i in range(len(id_pfrs)):
                    if i < len(id_pfrs)-1:
                        f.write('%d       input.pfr.%d.dic\n' % (id_pfrs[i], id_pfrs[i]))
                    else:
                        f.write('%d       input.pfr.%d.dic;\n' % (id_pfrs[i], id_pfrs[i]))
            # Write internal connections
            ncounts = 0  # Connections counter
            f.write('@InternalConnections     ')
            for i in range(Nr):
                for j in range(Nr):
                    if self.MassFlowrates[i,j] != 0 and ncounts < Nconns - 1:
                        f.write('%d     %d     %e \n' % (i, j, self.MassFlowrates[i,j]))
                        ncounts += 1
                    elif self.MassFlowrates[i,j] != 0 and ncounts == Nconns - 1:
                        f.write('%d     %d     %e;\n' % (i, j, self.MassFlowrates[i,j]))
                        ncounts += 1
            # Write inlet streams
            f.write('@InputStreams     ')
            for i in range(len(InletList)):
                if i < len(InletList) - 1:
                    f.write('%d     %e \n' % (InletList[i], self.Rlist[InletList[i]].Mf))
                elif i == len(InletList) - 1:
                    f.write('%d     %e; \n' % (InletList[i], self.Rlist[InletList[i]].Mf))
            # Write outlet streams
            f.write('@OutputStreams     ')
            for i in range(len(OutletList)):
                ri = OutletList[i]
                mi = np.abs(np.sum(self.MassFlowrates[ri,:]) - np.sum(self.MassFlowrates[:,ri]))
                if i < len(OutletList) - 1:
                    f.write('%d     %e \n' % (ri, mi))
                elif i == len(OutletList) - 1:
                    f.write('%d     %e; \n' % (ri, mi))
            # Write species to monitor
            f.write('@SpeciesToMonitor     ')
            for i in range(len(self.SpeciesMonitor)):
                f.write('%s\t' % self.SpeciesMonitor[i])
            f.write(';\n')
            # Verbosity Level
            f.write('@VerbosityLevel     1;\n')
            f.write('} \n \n')

            # Kinetic mechanism dictionary
            f.write('Dictionary kinetic-mechanism \n { \n')
            f.write('@Kinetics          %s;\n' % self.KinFile)
            f.write('@Thermodynamics    %s;\n' % self.ThermoFile)
            f.write('@Output     kinetics;\n')

            f.write('} \n')

        # Return simulation folder
        self.SimFolder = FolderName
        return self
        
    # Command to run simulation
    def RunSimulation(self, netsmoke_path, verbose=False):
        """
        Command to run simulation.

        Args:
            netsmoke_path (str): Path to the NetSMOKE++ executable.

        Raises:
            ValueError: If simulation folder is not specified as attribute in reactor network class.
        """

        # Check if SimFolder exists
        if os.path.exists(self.SimFolder) == False:
            raise ValueError('Simulation folder is not specified as attribute in reactor network class')
        
        # Create Run.sh file for simulation running
        runfile = self.SimFolder + '/' + 'Run.sh'
        with open(runfile, 'w') as f:
            f.write(netsmoke_path+'/SeReNetSMOKEpp.sh --input input.dic')
        
        # Create string for command to be executed
        os.chdir(self.SimFolder)
        command = 'sh Run.sh'
        log_file = "netlog.txt"
        if verbose == False:
            os.system(f"{command} >> {log_file}")
        else:
            os.system(command)
        # Return to working directory
        os.chdir(self.WorkingDir)

    def ExtractSingleOutput(self, Rid):

        """
        Extracts the output of a specific reactor and returns it as a Cantera Quantity object.

        Args:
            Rid (int): The index of the reactor.

        Returns:
            Cantera Quantity: The extracted reactor output.

        Raises:
            ValueError: If the simulation folder is not specified as an attribute in the reactor network class.
        """

        # Check if SimFolder exists
        if os.path.exists(self.SimFolder) == False:
            raise ValueError('Simulation folder is not specified as attribute in reactor network class')
        
        # Extract data for reactor Rid
        filename = self.SimFolder + '/Output/Reactor.' + str(Rid) + '/Output.out'
        df = pd.read_csv(filename, sep='\s+')

        # Export T, P
        T = df.values[-1,4]
        P = df.values[-1,5]

        # In columns, locate all the strings that have '_x" inside
        sp_list = []
        X_val   = []
        for i, col in enumerate(df.columns):
            ss = col.split('_x')
            if len(ss) > 1:
                sp_list.append(ss[0])
                X_val.append(df.values[-1,i])

        # Create solution
        gas = ct.Solution(self.CanteraMech)
        # Extract number of species
        ns = gas.n_species
        sp_list = gas.species_names
        X_old = gas.X
        X_new = X_old
        for i, sp in enumerate(sp_list):
            X_new[gas.species_index(sp)] = X_val[i]

        # Set new gas state
        gas.TPX = T, P, X_new
        # Create quantity object with mass from reactor
        M = ct.Quantity(gas, constant='HP', mass=self.Rlist[Rid].Mf)

        return M
    
    def ExtractOutputs(self):

        """
        Extracts the output of all reactors in the network and returns them as a list of Cantera Quantity objects.

        Returns:
            list: A list of Cantera Quantity objects representing the output of each reactor.

        Raises:
            ValueError: If the simulation folder is not specified as an attribute in the reactor network class.
        """

        # Check if SimFolder exists
        if os.path.exists(self.SimFolder) == False:
            raise ValueError('Simulation folder is not specified as attribute in reactor network class')
        
        # Initialize list of cantera quantities
        Mlist = []
        
        # Extract data for reactor Rid
        for i in range(self.Nr):
            M = self.ExtractSingleOutput(i)
            Mlist.append(M)
    
        return Mlist
    
    def ExtractSingleInput(self, Rid):

        '''This function will extract the input of the specific
        reactor Reactor.Rid that is specified and return it as a 
        cantera quantity object'''

        # Check if SimFolder exists
        if os.path.exists(self.SimFolder) == False:
            raise ValueError('Simulation folder is not specified as attribute in reactor network class')
        
        # Extract data for reactor Rid
        filename = self.SimFolder + '/Output/Reactor.' + str(Rid) + '/log.inlet'
        df = pd.read_csv(filename, sep='\s+')

        # Export T, P
        T = df.values[-1,1]

        # Extract data for reactor Rid
        filename2 = self.SimFolder + '/Output/Reactor.' + str(Rid) + '/Output.out'
        df2 = pd.read_csv(filename2, sep='\s+')
        P = df2.values[-1,5]

        # In columns, locate all the strings that have '_x" inside
        sp_list = []
        Y_val   = []
        for i, col in enumerate(df.columns[2:]):
            ss = col.split('(')
            if len(ss) > 1:
                sp_list.append(ss[0])
                Y_val.append(df.values[-1,i])

        # Create solution
        gas = ct.Solution(self.CanteraMech)
        # Extract number of species
        ns = gas.n_species
        sp_list = gas.species_names
        Y_old = gas.Y
        Y_new = Y_old
        for i, sp in enumerate(sp_list):
            Y_new[gas.species_index(sp)] = Y_val[i]

        # Set new gas state
        gas.TPY = T, P, Y_new
        # Create quantity object with mass from reactor
        M = ct.Quantity(gas, constant='HP', mass=self.Rlist[Rid].Mf)

        return M
    
    def ExtractInputs(self):

        """
        Extracts the input of the specific reactor specified by Rid and returns it as a Cantera Quantity object.

        Args:
            Rid (int): The ID of the reactor for which to extract the input.

        Returns:
            ct.Quantity: A Cantera Quantity object representing the input of the specified reactor.

        Raises:
            ValueError: If the simulation folder is not specified as an attribute in the reactor network class.
        """

        # Check if SimFolder exists
        if os.path.exists(self.SimFolder) == False:
            raise ValueError('Simulation folder is not specified as attribute in reactor network class')
        
        # Initialize list of cantera quantities
        Mlist = []
        
        # Extract data for reactor Rid
        for i in range(self.Nr):
            M = self.ExtractSingleInput(i)
            Mlist.append(M)
    
        return Mlist
    
    def ExtractOutputSingleVar(self, varname):

        """
        Extracts a single variable for each reactor and returns it as a NumPy array.

        Args:
            varname (str): The name of the variable to extract.

        Returns:
            numpy.ndarray: An array containing the values of the specified variable for each reactor.

        Raises:
            ValueError: If the simulation folder is not specified as an attribute in the reactor network class,
                        or if the variable is not found in the output file columns.
        """

        # Check if SimFolder exists
        if os.path.exists(self.SimFolder) == False:
            raise ValueError('Simulation folder is not specified as attribute in reactor network class')
        
        # Find the desired variable
        if varname == "T":
            varid = 4
        elif varname == "time" or varname == "tau":
            varid = 0
        else:
            f0 = self.SimFolder + '/Output/Reactor.0/Output.out'
            df = pd.read_csv(f0, sep='\s+')
            for j, col in enumerate(df.columns):
                ss = col.split('_x')
                if len(ss) > 0:
                    if varname == ss[0]:
                        varid = j
            if varid == 0:
                raise ValueError("The variable was not found in columns!")

        # Initialize output
        yout = np.zeros(self.Nr)
        for i in range(self.Nr):    
            filename = self.SimFolder + '/Output/Reactor.' + str(i) + '/Output.out'
            df = pd.read_csv(filename, sep='\s+')
            yout[i] = df.values[-1,varid]

        return yout
    
    def ExtractOutputMixtureFraction(self, fuel, oxidizer):

        """"
        Extracts the mixture fraction for each reactor and returns it as a list.

        Args:
            fuel (str): The name of the fuel species.
            oxidizer (str): The name of the oxidizer species.

        Returns:
            list: A list containing the mixture fraction for each reactor.

        Raises:
            ValueError: If the simulation folder is not specified as an attribute in the reactor network class.
        """

        # Extract all the outputs as a cantera quantity object
        Mlist = self.ExtractOutputs()
        Zlist = []
        for i in range(len(Mlist)):
            Zlist.append(Mlist[i].mixture_fraction(fuel=fuel, oxidizer=oxidizer, element="Bilger"))

        return Zlist

    def ExtractOutputEquivalenceRatio(self, fuel, oxidizer):

        """
        Extracts the equivalence ratio for each reactor and returns it as a list.

        Args:
            fuel (str): The name of the fuel species.
            oxidizer (str): The name of the oxidizer species.

        Returns:
            list: A list containing the equivalence ratio for each reactor.

        Raises:
            ValueError: If the simulation folder is not specified as an attribute in the reactor network class.
        """

        # Extract all the outputs as a cantera quantity object
        Mlist = self.ExtractOutputs()
        philist = []
        for i in range(len(Mlist)):
            philist.append(Mlist[i].equivalence_ratio(fuel=fuel, oxidizer=oxidizer))

        return philist

'''Class for generalized CRN'''
class GeneralCRN:

    def __init__(self, KinFile, Thermofile, CanteraMech="gri30.cti"):

        self.Kinfile_ = KinFile
        self.Thermofile_ = Thermofile
        self.CanteraMech_ = CanteraMech

    # Define function to initialize input layer
    def InputLayer(self, Inputs, taui=1e-6):

        # Number of inputs
        nin = len(Inputs)

        # Create each node as input node
        nodes = []
        locs  = []
        # Initialize list of reference mixtures and parameters
        ref_mixtures = []
        param_list = []
        for i in range(nin):
            nodes.append(("I"+str(i),  "PSR"))
            locs.append((1, i))
            # For inlets, the reference input is the Cantera quantity object
            ref_mixtures.append(Inputs[i])
            # Create the dictionary for standard inlet parameters
            params = {'tau':taui, 'Mf':Inputs[i].mass, "volume":None}
            param_list.append(params)

        # Create the attribute of layers
        self.layers_ = {"input":nodes}
        self.node_locations_ = locs
        # Initialize connections
        self.connections_ = []
        self.splitting_ratios_ = []
        # Initialze reference mixture and parameter list
        self.ref_mixtures_ = ref_mixtures
        self.param_list_ = param_list

        return self
    
    def AddConnectedLayer(self, nr, rtype="PSR", params=None, param_list=None, 
                          taui=1e-3, Li=1.0, Di=0.01, ref_mixture=None):

        '''This function will add a fully cinnected layer of reactors.
        The default reactor type is PSR, otherwise you can specify PFR via the
        keyword rtype. The keyword params is a dictionary containing parameters that
        will be used to initialize each reactor in the layer. Otherwise, if you want to
        pass specific parameters to each reactor, you can do it using the list param_list,
        which is a list of dictionaries, and each dictionary contains the parameters
        for each reactor'''

        # Check if params is None
        if params == None and rtype == "PSR" and param_list==None:
            pars = {'tau':taui, 'volume':None}
            param_list = []
            for i in range(nr):
                param_list.append(pars)
        elif params == None and rtype == "PFR" and param_list == None:
            pars = {'L':Li, 'D':Di, 'volume':None}
            param_list = []
            for i in range(nr):
                param_list.append(pars)

        # If param_list is given, none of the two if willl be executed and param_list will be used

        # Count existing nodes
        total_nodes = sum(len(nodes) for nodes in self.layers_.values())

        # Create new connections with all the nodes in the previous layer
        n_last = len(list(self.layers_.values())[-1])
        last_nodes = list(self.layers_.values())[-1]
        last_layer_name = list(self.layers_.keys())[-1]
        # Check if last_layer_name is input
        if last_layer_name == "input":
            ns = "H0"
            nl = 0
            layername = 'hidden0'
        else:
            nl = int(last_layer_name[6:]) + 1
            ns = "H" + str(nl)
            layername = 'hidden'+str(nl)

        newconns = []
        splitting_ratios = []
        split_ratio = round(1/nr, 4)
        for i in range(n_last):
            for j in range(nr):
                
                # Name of the new node
                newnode_name = ns + str(j)
                
                # Add new connection between previous layer and new layer
                newconns.append((last_nodes[i][0], newnode_name))
                splitting_ratios.append(split_ratio)

        # Create new hidden layer in the layers
        hidden_layer = []
        newlocs = []
        for i in range(nr):
            newnode_name = ns + str(i)
            hidden_layer.append((newnode_name, rtype))
            newlocs.append((int(nl+1)*2+1, i))

        # Update layers
        self.layers_[layername] = hidden_layer
        # Update nodes locations
        for item in newlocs:
            self.node_locations_.append(item)
        # Update nodes connections
        for item in newconns:
            self.connections_.append(item)
        # Update splitting rations
        for item in splitting_ratios:
            self.splitting_ratios_.append(item)
        # Update parameters list
        for item in param_list:
            self.param_list_.append(item)

        # Update reference mixtures
        if ref_mixture == None:
            gas = ct.Solution(self.CanteraMech_)
            gas.TPX = 300.0, 101325.0, "O2:0.21, N2:0.79"
            ref_mixture = gas
        for i in range(nr):
            self.ref_mixtures_.append(ref_mixture)
        
        return self
    
    def AddOutputLayer(self, nr, rtype="PSR", params=None, param_list=None, 
                          taui=1e-3, Li=1.0, Di=0.01, ref_mixture=None):
        
        # Check if params is None
        if params == None and rtype == "PSR" and param_list==None:
            pars = {'tau':taui, 'volume':None}
            param_list = []
            for i in range(nr):
                param_list.append(pars)
        elif params == None and rtype == "PFR" and param_list == None:
            pars = {'L':Li, 'D':Di, 'volume':None}
            param_list = []
            for i in range(nr):
                param_list.append(pars)

        # Count existing nodes
        total_nodes = sum(len(nodes) for nodes in self.layers_.values())

        # Create new connections with all the nodes in the previous layer
        n_last = len(list(self.layers_.values())[-1])
        last_nodes = list(self.layers_.values())[-1]
        last_layer_name = list(self.layers_.keys())[-1]
        nl = int(last_layer_name[6:]) + 1
        ns = "O"
        layername = 'output'

        newconns = []
        splitting_ratios = []
        split_ratio = round(1/nr, 4)
        for i in range(n_last):
            for j in range(nr):
                
                # Name of the new node
                newnode_name = ns + str(j)
                
                # Add new connection between previous layer and new layer
                newconns.append((last_nodes[i][0], newnode_name))
                splitting_ratios.append(split_ratio)

        # Create new hidden layer in the layers
        hidden_layer = []
        newlocs = []
        for i in range(nr):
            newnode_name = ns + str(i)
            hidden_layer.append((newnode_name, rtype))
            newlocs.append((int(nl+1)*2+1, i))

        # Update layers
        self.layers_[layername] = hidden_layer
        # Update nodes locations
        for item in newlocs:
            self.node_locations_.append(item)
        # Update nodes connections
        for item in newconns:
            self.connections_.append(item)
        # Update splitting rations
        for item in splitting_ratios:
            self.splitting_ratios_.append(item)
        # Update parameters list
        for item in param_list:
            self.param_list_.append(item)

        # Update reference mixtures
        if ref_mixture == None:
            gas = ct.Solution(self.CanteraMech_)
            gas.TPX = 300.0, 101325.0, "O2:0.21, N2:0.79"
            ref_mixture = gas
        for i in range(nr):
            self.ref_mixtures_.append(ref_mixture)

        return self
    
    def CreateGraph(self, plot=False):
        '''This function creates the graph of the CRN using the networkX library'''
        G = nx.Graph()
        # Add the nodes
        i=0
        for layer, nodes in self.layers_.items():
            for node, node_type in nodes:
                if node_type == "PFR":
                    shape="box"
                    G.add_node(node, type=node_type, shape=shape, ref_mixture=self.ref_mixtures_[i], parameters=self.param_list_[i], Rid=i)
                else:
                    G.add_node(node, type=node_type, ref_mixture=self.ref_mixtures_[i], parameters=self.param_list_[i], Rid=i)
                i+=1

        # Add the connections
        for i in range(len(self.connections_)):
            G.add_edge(self.connections_[i][0], self.connections_[i][1], weight=float(self.splitting_ratios_[i]))

        # Plot if true
        if plot == True:
            # Create dictionary of nodes positions
            locations = {}
            nodes     = list(G.nodes())
            for i in range(len(nodes)):
                locations[nodes[i]] = self.node_locations_[i]
            node_colors = {'PSR': 'lightblue', 'PFR': 'lightgreen'}
            node_type = nx.get_node_attributes(G, 'type')
            nx.draw(G, pos=locations, with_labels=True, node_color=[node_colors[node_type[node]] for node in G.nodes()],
                node_size=1500, edge_color='gray')
            plt.show()

        # Add graph to self
        self.G_ = G

        return self
    
    def GetMassFlowrates(self):

        # Get total number of reactors
        nr_tot = len(list(self.G_.nodes))
        self.nr_tot_ = nr_tot

        # Initialize the alpha matrix of splitting ratios
        alpha = np.zeros((nr_tot, nr_tot))

        # Assign split
        # Assign mass flow rates based on reactor indices
        for edge in self.G_.edges(data=True):
            node1, node2, weight = edge[0], edge[1], edge[2]['weight']
            Rid1 = self.G_.nodes[node1]['Rid']
            Rid2 = self.G_.nodes[node2]['Rid']
            alpha[Rid1, Rid2] = weight

        # Create array of boundary conditions
        bc = np.zeros(nr_tot)
        for i, node in enumerate(self.G_.nodes):
            if node.startswith('I'):
                mass = self.G_.nodes[node]['ref_mixture'].mass
                bc[i] = mass

        # Solve the linear system
        A = np.eye(nr_tot) - alpha.T
        M = np.linalg.solve(A,bc)

        mass_flowrates = np.zeros((nr_tot, nr_tot))
        for i in range(len(M)):
            mass_flowrates[i,:] = alpha[i,:] * M[i]

        # Get the outlet mass flowrates
        bc_out = np.zeros(self.nr_tot_)
        for i, node in enumerate(self.G_.nodes):
            if node.startswith('O'):
                bc_out[i] = sum(mass_flowrates[:,i]) - sum(mass_flowrates[i,:])
                if bc_out[i] < 0:
                    raise ValueError("Negative outlet mass flowrate detected in reactor %d" % i)
        

        self.mass_flowrates_ = mass_flowrates
        self.bc_in_ = bc
        self.bc_out_ = bc_out

        return self
    
    def CheckMassBalance(self):

        reactor_in  = np.zeros(self.nr_tot_)
        reactor_out = np.zeros(self.nr_tot_)
        local_imbalance = np.zeros(self.nr_tot_)

        # Mass balance of the network
        for j in range(self.nr_tot_):
            reactor_in[j]  = sum(self.mass_flowrates_[:,j])             # Mass inflow from other reactors
            reactor_out[j] = sum(self.mass_flowrates_[j,:])             # Mass outflow towards other reactors
            ext_input = self.bc_in_[j]                                  # Mass inflow from inlets
            ext_output = self.bc_out_[j]                               # Mass inflow from outlets
            inl = reactor_in[j] + ext_input                             # Total input
            out = reactor_out[j] + ext_output                           # Total output
            local_imbalance[j] = abs(inl-out)/max(inl,out)              # Local imbalance (relative %)

        # Total input and total output in the network
        total_in  = sum(reactor_in)
        total_out = sum(reactor_out)

        global_imbalance = abs(total_in - total_out)/max(total_in, total_out)

        # Check if the mass balance is globally satisfied
        mass_check = True
        if global_imbalance > 1e-12:
            raise ValueError('Check internal consistency of the network')
        
        else:
            # Check if the local mass balance is satisfied
            for j in range(self.nr_tot_):
                if local_imbalance[j] > 0.01: 
                    mass_check = False
                    stri = "Mass balance over 1 percent detected in reactor " + str(j)
                    raise ValueError(stri)

        if mass_check == True:
            print("Global mass imbalance and local mass imbalance within tolerances")
            print("Global imbalance = ", global_imbalance*100, ' percent')

        return self

    def WriteInputs(self):

        # Initialize reactor list
        rlist = []

        # Iterate over the nodes of the CRN to define the single reactors
        for i, node in enumerate(self.G_.nodes):
            rtype = self.G_.nodes[node]['type']
            parameters = self.G_.nodes[node]['parameters']
            mixture = self.G_.nodes[node]['ref_mixture']
            rid = self.G_.nodes[node]['Rid']

            # Get mass flowrate
            Mf = sum(self.mass_flowrates_[:,i])

            # Create reactor object
            if rtype == "PSR":
                isinput = False
                isoutput = False
                # Check if reactor is input or output
                if node.startswith('I'):
                    isinput=True
                    # Get mass flowrate
                    Mf = mixture.mass
                if node.startswith('O'):
                    isoutput=True 

                # Check if volume is a parameter of parameters
                if parameters["volume"] == None:
                    r = Reactor(Rtype=rtype, tau=parameters["tau"], Mf=Mf, CanteraMech=self.CanteraMech_, InletMixture=mixture,
                                InitialStatus=mixture, isinput=isinput, isoutput=isoutput)
                else:
                    r = Reactor(Rtype=rtype, volume=parameters["volume"], Mf=Mf, CanteraMech=self.CanteraMech_, InletMixture=mixture,
                                InitialStatus=mixture, isinput=isinput, isoutput=isoutput)
            # Create PFR reactor
            elif rtype == "PFR":

                isinput = False
                isoutput = False
                # Check if reactor is input or output
                if node.startswith('I'):
                    isinput=True
                if node.startswith('O'):
                    isoutput=True

                r = Reactor(Rtype=rtype, L=parameters["L"], D=parameters["D"], Mf=Mf, CanteraMech=self.CanteraMech_, InletMixture=mixture,
                                isinput=isinput, isoutput=isoutput)
                
            # Append reactor to list
            rlist.append(r)

        # Create reactor network object
        rn = ReactorNetwork(rlist, self.mass_flowrates_, KinFile=self.Kinfile_, ThermoFile=self.Thermofile_, CanteraMech=self.CanteraMech_)
        rn.WriteNetworkInput()

        # Add network to object
        self.rn_ = rn

        return self
    
    def RunCRN(self, netsmoke_path):

        self.rn_.RunSimulation(netsmoke_path)

        return self