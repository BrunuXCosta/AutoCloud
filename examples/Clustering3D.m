% AutoCloud v1.0 (previously referred as AutoClass)
% Developed by Bruno Sielly Jales Costa, Clauber Gomes Bezerra, Luiz
% Affonso Guedes and Plamen Angelov
% Commercial use not permitted
% Academic use only - with permission from authors
% Please cite the following papers:
% http://www.sciencedirect.com/science/article/pii/S0925231214013174
% http://ieeexplore.ieee.org/abstract/document/7502508/

clc;
clear;
addpath('../');

% Data set to be used
load('fisheriris.mat');
x = meas(:, 2:3);

% The only parameter required
m = 2;

% AutoCloud constructor
clusterer = AutoCloud('M', m, 'GenerateChart', true);

% Iteration over all data samples
for k = 1 : size(x, 1)
    [clusterer, ~, ~] = clusterer.addPoint(x(k, :));     
end