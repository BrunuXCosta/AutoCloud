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
X = load('../datasets/iris2.txt');

% 60% of the data is used for training, 40% for validation
size_train = 0.6 * size(X, 1);
Y = X(:, 5);
X = X(:, 1:4);
x_train = X(1:size_train, :);
x_test = X(size_train + 1:size(X, 1), :);
y_train = Y(1:size_train, :);
y_test = Y(size_train + 1:size(X, 1), :);

% The only parameter required
m = 3;

% AutoCloud constructor
classifier = AutoCloud('M', m, 'SimilarityMeasure', 'mahalanobis');

% Training
for k = 1 : size_train
    % All labels from training data are passed to the classifier
    [classifier, idx, ~] = classifier.addPoint(x_train(k, :), y_train(k));     
end

% Validating
output = zeros(size(X, 1) - size_train, 3);
target = zeros(size(X, 1) - size_train, 3);
for k = 1 : size(X, 1) - size_train    
    [classifier, idx, ~] = classifier.addPoint(x_test(k, :)); 
    output(k, idx) = 1;
    target(k, y_test(k)) = 1;
end

% Confusion Matrix
plotconfusion(target', output');