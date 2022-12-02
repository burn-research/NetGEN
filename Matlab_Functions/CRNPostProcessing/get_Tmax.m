function [Tmax] = get_Tmax(output_path)
% This function get the maximum temperature of a reactor network simulation

% Go into the output folder
cd(output_path);

Tmax = 0;

for i = 1: 100
    
    if isfolder(append('Reactor.', num2str(i))) == true
        cd(append('Reactor.', num2str(i)));
        data = importdata('Output.out');
        val = data.data;
        Tmax = max(Tmax, val(5));
        cd ../
    end
    
end

cd ../
        
    

end

