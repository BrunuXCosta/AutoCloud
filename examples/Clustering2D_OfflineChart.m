% AutoCloud v1.0 (previously referred as AutoClass)
% Developed by Bruno Sielly Jales Costa, Clauber Gomes Bezerra, Luiz
% Affonso Guedes and Plamen Angelov
% Commercial use not permitted
% Academic use only
% Please cite the following papers:
% http://www.sciencedirect.com/science/article/pii/S0925231214013174
% http://ieeexplore.ieee.org/abstract/document/7502508/

clc;
clear;
addpath('../');

% Data set to be used
x = load('../datasets/s2.txt');

% Max number of clouds to be plotted (only for plotting purposes)
maxClouds = 50;

dataSize = size(x, 1);
clusterer = AutoCloud('m', 2);
colorPick = hsv(maxClouds);
rng('default');
colorPick = colorPick(randperm(maxClouds), :);
colors = zeros(dataSize, 3);
for k = 1 : dataSize     
     % add new point (on-the-fly) to the classifier
     [clusterer, idx, ~] = clusterer.addPoint(x(k,:));     
     % get the cluster number with maximum membership     
     colors(k, :) = colorPick(idx, :);
     % get the center of the existing clusters so far     
     if (mod(k * 10, dataSize) == 0)
        fprintf('%d%%... ', floor(k*100/dataSize));
     end
end
centers = clusterer.getCenters();

scatter(x(:, 1), x(:, 2), 50, colors, '.');
hold on;
h = scatter(centers(:,1), centers(:,2), 100, 'k', 'o', 'LineWidth', 2);
drawnow;
hold off;