function [labels_new] = rewrite_fluent_labels(labels)
% This function transforms the fluent labels, e.g. '     mole-f ch4'
% in more readable form like 'CH4'

ns = length(labels);
labels_new = cell(1,ns);

for i = 1 : ns

    l = labels{i};
    ls = split(l);

    s = upper(extractAfter(ls, '-'));
    labels_new{i} = s{2};
end