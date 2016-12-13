classdef Cloud
    properties
        mu;
        var;
        n;
        covmat;
        name;
    end
    methods
        function obj = Cloud(varargin)
            obj.n = 0;
            obj.var = 0;
            obj.name = 'Default Class';
            obj.covmat = 0;
            if(nargin > 0)                
                obj.mu = varargin{1};
                obj.n = 1;        
                obj.covmat = zeros(size(obj.mu, 2), size(obj.mu, 2));
            end            
        end
        
        function obj = updateCloud(obj, mu, var, n, covmat)            
             obj.mu = mu;
             obj.var = var;
             obj.n = n;   
             obj.covmat = covmat;
        end        
        
        function obj = addPoint(obj, x)
            if (obj.n == 0)
                obj.n = 1;
                obj.mu = x;
                obj.var = 0;
                obj.covmat = zeros(size(obj.mu, 2), size(obj.mu, 2));
            else            
                obj.n = obj.n + 1;                
                obj.mu = ((obj.n - 1) / obj.n) * obj.mu + (1 / obj.n) * x;
                obj.var = ((obj.n - 1) / obj.n) * obj.var + (1 / (obj.n)) * norm(x - obj.mu) ^ 2;             
                obj.covmat = ((obj.n - 1) / (obj.n)) * obj.covmat + (1 / (obj.n)) * (x - obj.mu)' * (x - obj.mu);             
            end
        end
        
        function [zeta] = calculateZeta(obj, x, similarityMeasure)
            if (obj.n == 0)
                zeta = Inf;
                return
            end
            n_ = obj.n + 1;
            mu_ = ((n_ - 1) / n_) * obj.mu + (1 / n_) * x;
            var_ = ((n_ - 1) / n_) * obj.var + (1 / n_) * norm(x - mu_) ^ 2;            
            covmat_ = ((n_ - 1) / n_) * obj.covmat + (1 / n_) * (x - mu_)' * (x - mu_);
            if (strcmpi(similarityMeasure, 'euclidean')) %% Euclidean distance
                ksi = max((1 / n_) + ((mu_ - x) * (mu_ - x)') / (n_ * var_), 0.0001);
            else %% Mahalanobis distance
                ksi = max((1 / n_) + ((x - mu_) * pinv(covmat_) * (x - mu_)') / (n_ * size(mu_, 2)), 0.0001);
            end
            zeta = ksi / 2;           
        end
    end
end