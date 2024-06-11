function [split_ratio, mass_in, mass_out] = mass_split(mass_flowrates, bc_mass_flowrates)

[m,n] = size(mass_flowrates);

if m ~= n
    error('Matrix of mass flowrates is non squared');
end

% Get split ratio between the reactors
split_ratio = zeros(m, m);
for i = 1 : m
    if bc_mass_flowrates(i,2) == 0
        split_ratio(i,:) = mass_flowrates(i,:)/sum(mass_flowrates(i,:));
    elseif bc_mass_flowrates(i,2) ~= 0
        split_ratio(i,:) = mass_flowrates(i,:)/(sum(mass_flowrates(i,:)) - bc_mass_flowrates(i,2));
    end
end

% Get the inlet mass in each reactor
mass_in  = sum(mass_flowrates, 1)';
mass_out = sum(mass_flowrates, 2);

end

