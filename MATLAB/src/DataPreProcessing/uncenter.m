% This function provides a tool to center the data
% 
% INPUTS:
% 
% centered_data     = The matrix of centered data
% X_ave             = Matrix of previously substracted mean
%
% OUTPUT
%
% uncentered_data   = Matrix of uncentered data

function [uncentered_data] = uncenter(centered_data, X_ave)

uncentered_data = centered_data + X_ave;

%save rec_scaled_data uncentered_data