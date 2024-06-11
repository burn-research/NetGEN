function [output] = write_global_input(R, mass_flowrates, m_ext, opt)
% This function will create the global input file for netsmoke

% Write the global input files with the connections
file_id = fopen('input.dic', 'w');
fprintf(file_id, 'Dictionary ReactorNetwork \n { \n');
fprintf(file_id, '@KineticsPreProcessor      kinetic-mechanism; \n');
fprintf(file_id, '@MinIterations            5;   \n');
fprintf(file_id, '@MaxIterations            500;  \n');
fprintf(file_id, '@AtomicErrorThreshold     1e-4;\n');
fprintf(file_id, '@NonIsothermalErrorThreshold	1e-4; \n');
fprintf(file_id, '@MaxUnbalance     0.05; \n');

% Check how many PFR are present in R
pfr_count = 0;
for i = 1 : length(R)
    if strcmp(R{i}.Type, 'pfr') == true
        pfr_count = pfr_count + 1;
    end
end
fprintf('%d pfr are present \n \n', pfr_count);

% Check how many PSR are present in R
psr_count = 0;
for i = 1 : length(R)
    if strcmp(R{i}.Type, 'cstr') == true
        psr_count = psr_count + 1;
    end
end
fprintf('%d psr are present \n \n', psr_count);

% If there are PFR:
pfr_count_2 = pfr_count-1;
if pfr_count ~= 0
    fprintf(file_id, '@PlugFlowReactors     \n');
    
    for i = 1 : length(R)
        if strcmp(R{i}.Type, 'pfr') == true && pfr_count_2 ~= 0
            react_dic = append('input.pfr.', num2str(R{i}.id), '.dic');
            fprintf(file_id, '%d        %s  \n', R{i}.id, react_dic);
            pfr_count_2 = pfr_count_2 - 1;
        elseif strcmp(R{i}.Type, 'pfr') == true && pfr_count_2 == 0
            react_dic = append('input.pfr.', num2str(R{i}.id), '.dic');
            fprintf(file_id, '%d        %s;  \n', R{i}.id, react_dic);
        end
    end
end

% If there are PSR:
psr_count_2 = psr_count-1;
if psr_count ~= 0
    fprintf(file_id, '@PerfectlyStirredReactors     \n');
    
    for i = 1 : length(R)
        if strcmp(R{i}.Type, 'cstr') == true && psr_count_2 ~= 0
            react_dic = append('input.cstr.', num2str(R{i}.id), '.dic');
            fprintf(file_id, '%d        %s  \n', R{i}.id, react_dic);
            psr_count_2 = psr_count_2 - 1;
        elseif strcmp(R{i}.Type, 'cstr') == true && psr_count_2 == 0
            react_dic = append('input.cstr.', num2str(R{i}.id), '.dic');
            fprintf(file_id, '%d        %s;  \n', R{i}.id, react_dic);
        end
    end
end

% Write internal connections
fprintf(file_id, '@InternalConnections   ');
n_conn = length(find(mass_flowrates ~= 0));
count_conn = 1;
k = length(R);
for i = 1 : k
    for j = 1 : k
        if mass_flowrates(i,j) ~= 0
            if count_conn ~= n_conn
                fprintf(file_id, '                  %d   %d   %d \n', i-1, j-1, mass_flowrates(i,j));
                count_conn = count_conn + 1;
            elseif count_conn == n_conn
                fprintf(file_id, '                  %d   %d   %d; \n', i-1, j-1, mass_flowrates(i,j));
            end
        end
    end
end

% Write the inlet streams
fprintf(file_id, '@InputStreams   ');
n_in = length(find(m_ext(:,1) > 0));
in_count = 1;
for j = 1 : k
    if m_ext(j,1) ~= 0 && in_count ~= n_in
        fprintf(file_id, '%d %d \n', j-1,m_ext(j,1));
        in_count = in_count + 1;
    elseif m_ext(j,1) > 0 && in_count == n_in
        fprintf(file_id, '%d %d; \n', j-1, m_ext(j,1));
    end
end

% Write the outlet streams
fprintf(file_id, '@OutputStreams    ');
n_out = length(find(m_ext(:,2) > 0));   % Number of outlet reactors
out_count = 1;
for j = 1 : k
    if m_ext(j,2) > 0 && out_count ~= n_out
        fprintf(file_id, '%d %d \n', j-1, m_ext(j,2));
        out_count = out_count + 1;
    elseif m_ext(j,2) > 0 && out_count == n_out
        fprintf(file_id, '%d %d; \n', j-1, m_ext(j,2));
    end
end

% Print main options
if isfield(opt, 'SpeciesMonitor') == true
    sp = opt.SpeciesMonitor;
    st = '';
    for i = 1 : length(sp)
        st = append(st, '  ', sp{i});
    end
    fprintf(file_id, '@SpeciesToMonitor %s ; \n', st);
    fprintf(file_id, '@VerbosityLevel   1; \n } \n');
else
    fprintf(file_id, '@SpeciesToMonitor   CH4 H2 O2 CO2 CO NO NO2; \n');
    fprintf(file_id, '@VerbosityLevel   1; \n } \n');
end

% Print the kinetic mechanism dictionary
% Check if kinetic mechanism is specified
kin_mech = 'gri30';
if isfield(opt, 'KineticPath') == true
    fprintf(file_id, 'Dictionary kinetic-mechanism \n { \n');
    fprintf(file_id, '@Kinetics %s; \n', opt.Chemfile);
    fprintf(file_id,  '@Thermodynamics %s; \n', opt.Thermofile);
    fprintf(file_id, '@Output kinetics; \n } \n');
else
    switch kin_mech
        case 'gri30'
            fprintf('Kinetics folder not specified, GRI 3.0 will be used by default \n');
            fprintf(file_id, 'Dictionary kinetic-mechanism \n { \n @Kinetics ');
            fprintf(file_id, '../kinetics/GRI30/GRI30.CKI; \n');
            fprintf(file_id, '@Thermodynamics ');
            fprintf(file_id, '../kinetics/GRI30/thermo30.CKT; \n');
            fprintf(file_id, '@Output kinetics; \n }');
    
        case 'GAL'
            fprintf(file_id, 'Dictionary kinetic-mechanism \n { \n @Kinetics ');
            fprintf(file_id, '../kinetics/GAL/KIN.MECH; \n');
            fprintf(file_id, '@Thermodynamics ');
            fprintf(file_id, '../kinetics/GAL/THERMO.THERM; \n');
            fprintf(file_id, '@Output kinetics; \n }');
    end
end
    
mess = sprintf('File in format .dic available in the folder Input_files ready for NetSMOKE input');
disp(mess);
disp('Operation terminated');

output = true;

end

