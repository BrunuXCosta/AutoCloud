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
x = meas(:, 3:4);

% The only parameter required
m = 2;

% AutoCloud constructor
clusterer = AutoCloud('M', m);

colors = zeros(size(x, 1), 3);
% Iteration over all data samples
for k = 1 : size(x, 1)
    [clusterer, ~, membership] = clusterer.addPoint(x(k, :)); 
    % The color of each point is defined by its fuzzy membership to each
    % existing class
    color = membership;
    while (size(color, 1) < 3)
        color = [color; 0];
    end
    colors(k, :) = color;
end
% Colors normalized
colors = colors - repmat(min(colors), size(colors, 1), 1);
colors = colors ./ repmat(max(colors), size(colors, 1), 1);
centers = clusterer.getCenters();
scatter(x(:, 1), x(:, 2), 50, colors, '*');
hold on;
h = scatter(centers(:,1), centers(:,2), 100, 'k', 'o', 'LineWidth', 2);
drawnow;
hold off;