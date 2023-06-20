function [output] = write_reactor_input(R)
% This function will create a NetSMOKE input dictionary given the input
% reactor R which comprise several attributes

if isfield(R, 'Type') == false
    warning('Reactor type not specified, will be set to PSR');
    R.Type = 'cstr';
end

% Check the type of the reactor, if it is a cstr:
if strcmp(R.Type, 'cstr') == true
    if isfield(R, 'id') == true
        file_name = append('input.cstr.', num2str(R.id), '.dic');
    else
        file_name = 'input.dic';
    end

    file_id = fopen(file_name, 'w');
    fprintf(file_id, 'Dictionary PerfectlyStirredReactor \n');
    fprintf(file_id, '{ \n');

    % Kinetics pre-processor
    if isfield(R, 'id') == true
        fprintf(file_id, '      @KineticsFolder     dummy; \n');
    else 
        fprintf(file_id, '      @KineticsPreProcessor     kinetics; \n');
    end
    
    % Inlet status
    fprintf(file_id, '      @InletStatus        Inlet-Mixture; \n');

    % Initialize with initial status
    if isfield(R, 'InitialStatus')
        fprintf(file_id, '       @InitialStatus     Initial-Status; \n');
    end

    % Check option for ode settings
    if isfield(R, 'odesettings') == true
        fprintf(file_id, '       @OdeParameters       ode-settings; \n');
    end

    % Check option for sensitivity analysis
    if isfield(R, 'sensitivity') == true
        fprintf(file_id, '        @SensitivityAnalysis      sensitivity; \n');
    end

    % Check option for ROPA
    if isfield(R, 'ROPA') == true
        fprintf(file_id, '        @OnTheFlyROPA              ropa; \n');
    end
    
    % Check if the reactor is isothermal
    if isfield(R, 'isothermal') == true
        fprintf(file_id, '      @Type               Isothermal-ConstantPressure; \n');
        
    else
        % Non-isothermal options
        fprintf(file_id, '      @Type               NonIsothermal-ConstantPressure; \n');
        if isfield(R, 'Qdot') == true
            R.U = R.Qdot/(R.A*(R.Tout-R.Tenv));
        end
        if isfield(R, 'U') == true
            fprintf(file_id, '      @GlobalThermalExchangeCoefficient   %d W/m2/K; \n', R.U);
            fprintf(file_id, '      @ExchangeArea               %d m2; \n', R.A);
            fprintf(file_id, '      @EnvironmentTemperature     %d K;  \n', R.Tenv);
        end
    end

    % Check for kinetic corrections
    if isfield(R, 'KineticCorrections') == true
        fprintf(file_id, '      @KineticCorrections     kinetic-corrections; \n');
    end

    % Turbulent diffusivity
    if isfield(R, 'TurbDiff') == true
        fprintf(file_id, '      @TurbulentDiffusivity     %d kg/m/s; \n', R.TurbDiff);
    end

    % Check constraints on V, Mf or tau
    fcount = 0;
    if isfield(R, 'V') == true
        fprintf(file_id, '      @Volume             %d cm3; \n', R.V);
        fcount = fcount + 1;
    end
    if isfield(R, 'Mf') == true
        fprintf(file_id, '      @MassFlowRate      %d kg/s;   \n', R.Mf);
        fcount = fcount + 1;
    end
    if isfield(R, 'tau') == true
        fprintf(file_id, '      @ResidenceTime      %d s;   \n', R.tau);
        fcount = fcount + 1;
    end

    % Check if reactor attributes are under or overspecified
    if fcount == 0 || fcount == 1
        error('Specify at least two between mass flowrate, volume and tau...');
    elseif fcount == 3
        warning('Reactor is overspecified, remove one between mass flowrate, tau and volume...');
    end
    fprintf(file_id, '} \n');
    
    % Inlet mixture
    fprintf(file_id, '\n Dictionary Inlet-Mixture \n { \n');
    fprintf(file_id, '      @Pressure           %d Pa; \n', R.P);
    fprintf(file_id, '      @Temperature        %d K;  \n', R.T);
    if isfield(R, 'basis') == false
        fprintf(file_id, '      @Masses       ');
        for i = 1 : length(R.Y)
            if i ~= length(R.Y)
                fprintf(file_id, '          %s %s\n', extractBefore(R.Y{i}, ':'), extractAfter(R.Y{i}, ':'));
            else 
                fprintf(file_id, '          %s %s; \n', extractBefore(R.Y{i}, ':'), extractAfter(R.Y{i}, ':'));
                fprintf(file_id, '} \n');
            end
        end
    else
        if strcmp(R.basis, 'mol') || strcmp(R.basis, 'mole') == true
            fprintf(file_id, '      @Moles       ');
            for i = 1 : length(R.Y)
                if i ~= length(R.Y)
                    fprintf(file_id, '          %s %s\n', extractBefore(R.Y{i}, ':'), extractAfter(R.Y{i}, ':'));
                else 
                    fprintf(file_id, '          %s %s; \n', extractBefore(R.Y{i}, ':'), extractAfter(R.Y{i}, ':'));
                    fprintf(file_id, '} \n');
                end
            end
        elseif strcmp(R.basis, 'mass') || strcmp(R.basis, 'masses') == true
            fprintf(file_id, '      @Masses       ');
            for i = 1 : length(R.Y)
                if i ~= length(R.Y)
                    fprintf(file_id, '          %s %s\n', extractBefore(R.Y{i}, ':'), extractAfter(R.Y{i}, ':'));
                else 
                    fprintf(file_id, '          %s %s; \n', extractBefore(R.Y{i}, ':'), extractAfter(R.Y{i}, ':'));
                    fprintf(file_id, '} \n');
                end
            end
        end
    end

    % Initial status
    if isfield(R, 'InitialStatus')
        fprintf(file_id, '\n Dictionary Initial-Status \n { \n');
        fprintf(file_id, ' @Temperature 2000 K; \n');
        fprintf(file_id, ' @Pressure %d Pa; \n', R.P);
        fprintf(file_id, ' @Moles O2 0.21 \n N2 0.79; \n } \n');
    end
            
    % ODE settings
    if isfield(R, 'odesettings') == true
        fprintf(file_id, '\n Dictionary ode-settings \n { \n');
        fprintf(file_id, '@AbsoluteTolerance     %d; \n', R.abstol);
        fprintf(file_id, '@RelativeTolerance     %d; \n', R.reltol);
        fprintf(file_id, '}');
    end
    
    % Sensitivity Analysis
    if isfield(R, 'sensitivity') == true
        fprintf(file_id, '\n Dictionary sensitivity \n { \n');
        fprintf(file_id, '     @Type     kinetic-constants; \n');
        fprintf(file_id, '     @DenseSolver      Eigen; \n');
        fprintf(file_id, '     @DenseFullPivoting     false; \n');
        fprintf(file_id, '     @SubSteps     5; \n');
        if isfield(R, 'sensSpecies') == false
            fprintf(file_id, '     @Species NO; \n');
        else
            ss = '';
            for l = 1 : length(R.sensSpecies)
                if l ~= length(R.sensSpecies)
                    append(ss, ' ', R.sensSpecies{l}, ' ');
                else
                    append(ss, ' ', R.sensSpecies{l}, '; \n');
                end
            end
            disp(ss); pause;
            fprintf(file_id, append('@Species ', ss));
        end
        fprintf(file_id, '} \n;');
    end

    % Ropa analysis
    if isfield(R, 'ROPA') == true
        fprintf(file_id, '\n Dictionary ropa \n { \n');
        % Check for species
        if isfield(R, 'ropaSpecies') == false
            fprintf(file_id, '     @Species NO; \n');
        else
            ss = '';
            for l = 1 : length(R.ropaSpecies)
                if l ~= length(R.ropaSpecies)
                    append(ss, ' ', R.ropaSpecies{l}, ' ');
                else
                    append(ss, ' ', R.ropaSpecies{l}, '; \n');
                end
            end
            disp(ss); pause;
            fprintf(file_id, append('@Species ', ss));
        end
        % Check for reference species
        if isfield(R, 'ropaReference') == false
            fprintf(file_id, '     @ReferenceSpecies NO; \n');
        else
            ss = '';
            for l = 1 : length(R.ropaReference)
                if l ~= length(R.ropaReference)
                    append(ss, ' ', R.ropaReference{l}, ' ');
                else
                    append(ss, ' ', R.ropaReference{l}, '; \n');
                end
            end
            disp(ss); pause;
            fprintf(file_id, append('@ReferenceSpecies ', ss));
        end
        % Check for threshold
        if isfield(R, 'ropaThreshold') == true
            fprintf(file_id, '@Threshold     1e-3; \n');
        else
            fprintf(file_id, '@Threshold     %d; \n', R.ropaThreshold);
        end
        % Merge forward - backward
        if isfield(R, 'ropaMerge') == true
            fprintf(file_id, '@MergeForwardAndBackwardReactions     true; \n');
        end
        
        fprintf(file_id, '} \n');
    end




    % Kinetic dictionary
    if isfield(R, 'kinetic') == true
        fprintf(file_id, '\n Dictionary kinetics \n { \n');
        fprintf(file_id, append('@Kinetics ', R.mech, '; \n'));
        fprintf(file_id, append('@Thermodynamics ', R.therm, '; \n'));
        fprintf(file_id, '@Output Kinetics; \n');
        fprintf(file_id, '} \n');
    end

    % Kinetic corrections dictionary
    if isfield(R, 'KineticCorrections') == true
        fprintf(file_id, '\n Dictionary kinetic-corrections \n { \n');
        fprintf(file_id, '     @AverageTemperature     %d; \n', R.Tmean);
        fprintf(file_id, '     @TemperatureVariance    %d; \n', R.Tvar);
        fprintf(file_id, '     @Type beta; \n } \n');
    end
    
 
% If it is a PFR reactor:
else
    file_name = append('input.pfr.', num2str(R.id), '.dic');
    file_id = fopen(file_name, 'w');
    fprintf(file_id, 'Dictionary PlugFlowReactor \n');
    fprintf(file_id, '{ \n');
    fprintf(file_id, '      @KineticsFolder      dummy;        \n');

    % Check if it is isothermal
    if isfield(R, 'Isothermal')
        fprintf(file_id, '      @Type               Isothermal; \n');
    else
        fprintf(file_id, '      @Type               NonIsothermal; \n');
    end

    fprintf(file_id, '      @InletStatus        Inlet-Mixture; \n');
    fprintf(file_id, '      @ConstantPressure   true;          \n');

    % Check which properties are given
    if isfield(R, 'D') == true
        fprintf(file_id, '      @Diameter           %d mm;         \n', R.D);
    end
    if isfield(R, 'V')
        fprintf(file_id, '      @Volume             %d cm3;         \n', R.V);
    end
    if isfield(R, 'L') == true
        fprintf(file_id, '      @Length             %d mm;         \n', R.L);
    end
    if isfield(R, 'Mf') == true
        fprintf(file_id, '      @MassFlowRate       %d kg/s;       \n', R.Mf);
    end
    if isfield(R, 'Tau') == true
        fprintf(file_id, '      @ResidenceTime      %d s;          \n', R.Tau);
    end
    if isfield(R, 'U')
        fprintf(file_id, '      @GlobalThermalExchangeCoefficient      %d W/m2/K;          \n', R.U);
    end
    if isfield(R, 'Tenv')
        fprintf(file_id, '      @EnvironmentTemperature      %d K;          \n', R.Tenv);
    end
    if isfield(R, 'CrossOverPerimeter') 
        fprintf(file_id, '      @CrossSectionOverPerimeter      %d mm;          \n', R.D/4);
    end


    fprintf(file_id, '}                                        \n');

    % Inlet mixture
    fprintf(file_id, '\n Dictionary Inlet-Mixture \n { \n');
    fprintf(file_id, '      @Pressure           %d Pa; \n', R.P);
    fprintf(file_id, '      @Temperature        %d K;  \n', R.T);
    fprintf(file_id, '      @Masses       ');
    for i = 1 : length(R.Y)
        if i ~= length(R.Y)
            fprintf(file_id, '          %s %s\n', extractBefore(R.Y{i}, ':'), extractAfter(R.Y{i}, ':'));
        else 
            fprintf(file_id, '          %s %s; \n', extractBefore(R.Y{i}, ':'), extractAfter(R.Y{i}, ':'));
            fprintf(file_id, '} \n');
        end
    end
end

output = true;
end

