function [output] = write_netsmoke_inputs(mass_flowrates, bc_mass_flowrates, inlets, case_info, Vpsr, folder_name, Rtype, Lpfr, Dpfr)
% This function, given the inputs, will write the NetSMOKE input
% dictionaries associated with the conditions

% Number of reactors
[k, ~] = size(mass_flowrates);

% Calculate the thermal exchange parameters
Qdot = case_info.Qdot;      % W
Tout = case_info.Tout;      % K
A = 1;                      % m2
Tenv = 300;                 % K
U = Qdot/(A*(Tout - Tenv)); % W/m2 K

% Apply a correction factor for numerical stability of the solver
a = 10000;

% Calculate the mass flowrate that enters in each reactor
mass_in = zeros(k,1);
for j = 1 : k
    if bc_mass_flowrates(j) > 0
        mass_in(j) = sum(mass_flowrates(:,j)) + bc_mass_flowrates(j);
    else
        mass_in(j) = sum(mass_flowrates(:,j));
    end
end

% Create a folder
mkdir(folder_name)
cd(folder_name)

% Write the single reactor.dic files
for j = 1 : k
    
    % Rtype(j) = 0 is a non-isothermal adiabatic PSR
    if Rtype(j) == 0
        file_name = append('input.cstr.', num2str(j-1), '.dic');
        file_id = fopen(file_name, 'w');
        fprintf(file_id, 'Dictionary PerfectlyStirredReactor \n');
        fprintf(file_id, '{ \n');
        fprintf(file_id, '      @KineticsFolder      dummy; \n');
        fprintf(file_id, '      @Type               NonIsothermal-ConstantPressure; \n');
        fprintf(file_id, '      @InletStatus        Inlet-Mixture; \n');
        fprintf(file_id, '      @Volume             %d cm3; \n', Vpsr(j));
        fprintf(file_id, '      @MassFlowRate      %d kg/s;   \n', mass_in(j));
        fprintf(file_id, '} \n');
        fprintf(file_id, '\n Dictionary Inlet-Mixture \n { \n');
        fprintf(file_id, '      @Pressure           101325.0 Pa; \n');
        
        % If the reactor is an inlet reactor
        if bc_mass_flowrates(j,1) ~= 0
            fprintf(file_id, '      @Temperature        %d K;   \n', inlets{j}.T);
            fprintf(file_id, '      @MassFractions       ');
            for i = 1 : length(inlets{j}.Y)
                if i ~= length(inlets{j}.Y)
                    fprintf(file_id, '          %s %s\n', extractBefore(inlets{j}.Y{i}, ':'), extractAfter(inlets{j}.Y{i}, ':'));
                else 
                    fprintf(file_id, '          %s %s; \n', extractBefore(inlets{j}.Y{i}, ':'), extractAfter(inlets{j}.Y{i}, ':'));
                    fprintf(file_id, '} \n');
                end
            end
            
        % If it is not an inlet reactor    
        else
            fprintf(file_id, '      @Temperature        %d K;   \n', 1500);
            fprintf(file_id, '      @MassFractions      %d O2   \n', 0.233);
            fprintf(file_id, '                          %d N2;  \n', 0.767);
        end
        
        fprintf(file_id, '} \n');
                
    % Rtype(j) = 2 means that is a non-isothermal PSR with heat exchange
    elseif Rtype(j) == 2
        file_name = append('input.cstr.', num2str(j-1), '.dic');
        file_id = fopen(file_name, 'w');
        fprintf(file_id, 'Dictionary PerfectlyStirredReactor \n');
        fprintf(file_id, '{ \n');
        fprintf(file_id, '      @KineticsFolder      dummy; \n');
        fprintf(file_id, '      @Type               NonIsothermal-ConstantPressure; \n');
        fprintf(file_id, '      @GlobalThermalExchangeCoefficient   %d W/m2/K; \n', U);
        fprintf(file_id, '      @ExchangeArea               %d m2; \n', A);
        fprintf(file_id, '      @EnvironmentTemperature     %d K;  \n', Tenv);
        fprintf(file_id, '      @InletStatus        Inlet-Mixture; \n');
        fprintf(file_id, '      @Volume             %d cm3; \n', Vpsr(j));
        fprintf(file_id, '      @MassFlowRate      %d kg/s;   \n', mass_in(j));
        fprintf(file_id, '} \n');
        fprintf(file_id, '\n Dictionary Inlet-Mixture \n { \n');
        fprintf(file_id, '      @Temperature           %d K;   \n', 1500);
        fprintf(file_id, '      @Pressure       101325.0 Pa;   \n');
        fprintf(file_id, '      @MassFractions      O2 0.232 \n');
        fprintf(file_id, '                          N2 0.768;\n');
        fprintf(file_id, '} \n');
    
    % Rtype(j) = 1 means that the reactor is a PFR    
    elseif Rtype(j) == 1
        file_name = append('input.pfr.', num2str(j-1), '.dic');
        file_id = fopen(file_name, 'w');
        fprintf(file_id, 'Dictionary PlugFlowReactor \n');
        fprintf(file_id, '{ \n');
        fprintf(file_id, '      @KineticsFolder      dummy;        \n');
        fprintf(file_id, '      @Type               NonIsothermal; \n');
        fprintf(file_id, '      @InletStatus        Inlet-Mixture; \n');
        fprintf(file_id, '      @ConstantPressure   true;          \n');    
        fprintf(file_id, '      @Diameter           %d mm;         \n', Dpfr(j));
        fprintf(file_id, '      @Length             %d mm;         \n', Lpfr(j));
        fprintf(file_id, '      @MassFlowRate       %d kg/s;       \n', mass_in(j));
        fprintf(file_id, '}                                        \n');
        fprintf(file_id, '\n Dictionary Inlet-Mixture \n { \n');
        fprintf(file_id, '      @Temperature           %d K;   \n', 1500);
        fprintf(file_id, '      @Pressure       101325.0 Pa;   \n');
        fprintf(file_id, '      @MassFractions      O2 0.232 \n');
        fprintf(file_id, '                          N2 0.768;\n');
        fprintf(file_id, '} \n');                
    end
end

% Write the global input files with the connections
file_id = fopen('input.dic', 'w');
fprintf(file_id, 'Dictionary ReactorNetwork \n { \n');
fprintf(file_id, '@KineticsPreProcessor      kinetic-mechanism; \n');
fprintf(file_id, '@MinIterations            5;   \n');
fprintf(file_id, '@MaxIterations            500;  \n');
fprintf(file_id, '@AtomicErrorThreshold     1e-3;\n');
fprintf(file_id, '@NonIsothermalErrorThreshold	1e-3; \n');
fprintf(file_id, '@MaxUnbalance     0.01; \n');

pfr = find(Rtype == 1);
psr = find(Rtype == 0 | Rtype == 2);

fprintf(file_id, '@PlugFlowReactors     \n');
for j = 1 : length(pfr)
    if j ~= length(pfr)
            react_dic = append('input.pfr.', num2str(pfr(j)-1), '.dic');
            fprintf(file_id, '%d        %s  \n', pfr(j)-1, react_dic);
    else
            react_dic = append('input.pfr.', num2str(pfr(j)-1), '.dic');
            fprintf(file_id, '%d        %s;  \n', pfr(j)-1, react_dic);
    end
end

fprintf(file_id, '@PerfectlyStirredReactors     \n');
for j = 1 : length(psr)
    if j ~= length(psr)
        react_dic = append('input.cstr.', num2str(psr(j)-1), '.dic');
        fprintf(file_id, '%d        %s  \n', psr(j)-1, react_dic);
    else
        react_dic = append('input.cstr.', num2str(psr(j)-1), '.dic');
        fprintf(file_id, '%d        %s;  \n', psr(j)-1, react_dic);
    end 
end

% Remember to apply a correction factor to have connections flowrate of the
% order of unity
fprintf(file_id, '@InternalConnections   ');
n_conn = length(find(mass_flowrates ~= 0));
count_conn = 1;
for i = 1 : k
    for j = 1 : k
        if mass_flowrates(i,j) ~= 0
            if count_conn ~= n_conn
                fprintf(file_id, '                  %d   %d   %d \n', i-1, j-1, a*mass_flowrates(i,j));
                count_conn = count_conn + 1;
            elseif count_conn == n_conn
                fprintf(file_id, '                  %d   %d   %d; \n', i-1, j-1, a*mass_flowrates(i,j));
            end
        end
    end
end
    
% Print the inlet streams
fprintf(file_id, '@InputStreams   ');
n_in = length(find(bc_mass_flowrates(:,1) > 0));
in_count = 1;
for j = 1 : k
    if bc_mass_flowrates(j,1) ~= 0 && in_count ~= n_in
        fprintf(file_id, '%d %d \n', j-1, bc_mass_flowrates(j,1)*a);
        in_count = in_count + 1;
    elseif bc_mass_flowrates(j,1) > 0 && in_count == n_in
        fprintf(file_id, '%d %d; \n', j-1, bc_mass_flowrates(j,1)*a);
    end
end

% Print the outlets
fprintf(file_id, '@OutputStreams    ');
n_out = length(find(bc_mass_flowrates(:,2) < 0));   % Number of outlet reactors
out_count = 1;
for j = 1 : k
    if bc_mass_flowrates(j,2) < 0 && out_count ~= n_out
        fprintf(file_id, '%d %d \n', j-1, -bc_mass_flowrates(j,2)*a);
        out_count = out_count + 1;
    elseif bc_mass_flowrates(j,2) < 0 && out_count == n_out
        fprintf(file_id, '%d %d; \n', j-1, -bc_mass_flowrates(j,2)*a);
    end
end

% Print main options
fprintf(file_id, '@SpeciesToMonitor NO NO2 NH3; \n');
fprintf(file_id, '@VerbosityLevel   1; \n } \n');

% Print the kinetic mechanism dictionary
fprintf(file_id, 'Dictionary kinetic-mechanism \n { \n @Kinetics ');
fprintf(file_id, '../kinetics/Ammonia/NH3_chem.CKI; \n');
fprintf(file_id, '@Thermodynamics ');
fprintf(file_id, '../kinetics/Ammonia/NH3_thermo.CKT; \n');
fprintf(file_id, '@Output kinetics; \n }');

cd ../
    
mess = sprintf('File in format .dic available in the folder Input_files ready for NetSMOKE input');
disp(mess);
disp('Operation terminated');

output = true;


end

