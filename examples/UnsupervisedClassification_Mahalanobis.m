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
x = load('../datasets/iris.txt');
y = x(:, 5);
x = x(:, 1:4);

% The only parameter required
m = 3;

% AutoCloud constructor
classifier = AutoCloud('M', m, 'SimilarityMeasure', 'mahalanobis');

% Iteration over all data samples
output = zeros(size(x, 1), 3);
target = zeros(size(x, 1), 3);
for k = 1 : size(x, 1)
    % No labels are known, no training is performed
    [classifier, idx, ~] = classifier.addPoint(x(k, :)); 
    output(k, idx) = 1;
    target(k, y(k)) = 1;
end

% Confusion Matrix
plotconfusion(target', output');