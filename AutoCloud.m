classdef AutoCloud
    properties(Access = public)
        k;
        cloudList;                       
    end
    properties(Access = private)
        m;
        intersectionList;
        intersectionMatrix;
        contMerge;
        dimension;
        membershipList;
        echo;
        similarityMeasure;
        generateChart;
        generateVideo;
        autoMerge;
        lHandle;
        cHandle;
        delay;
        colorPick;
        colors;
    end
    methods
        function obj = AutoCloud(varargin)
            if ~(nargin == 1 || (nargin >=2 && mod(nargin, 2) == 0))
               error('Wrong parameters. AutoCloud could not be initialized.'); 
            end            
            obj.k = 0;
            obj.m = 1;
            c = Cloud();
            c.name = 'Class 1';            
            obj.cloudList = c;
            obj.intersectionList = 0;
            obj.intersectionMatrix = 0;
            obj.contMerge = 0;
            obj.dimension = 0; 
            obj.similarityMeasure = 'Euclidean';
            obj.echo = false;
            obj.generateChart = false;
            obj.generateVideo = false;
            obj.autoMerge = true;
            obj.delay = 0.01;            
            obj.colorPick = hsv(256);
            %rng('default');
            %obj.colorPick = obj.colorPick(randperm(50), :);
            obj.colors = [];
            if (nargin == 1)
                obj.m = varargin{1};
            else
                for i = 1 : 2 : nargin
                    parameter = varargin{i};
                    if (~ischar(parameter))
                        error(strcat('Invalid parameter name (parameter #', num2str(i), ').'));
                    end
                    if (strcmpi(parameter, 'M'))
                        if (~isfloat(varargin{i + 1}) || varargin{i + 1} <= 0)
                            error(strcat('Invalid parameter value (parameter #', num2str(i + 1), ').'));
                        end
                        obj.m = varargin{i + 1};                    
                    elseif (strcmpi(parameter, 'SimilarityMeasure'))
                        if (~ischar(varargin{i + 1}) || ~strcmpi(varargin{i + 1}, 'euclidean') && ~strcmpi(varargin{i + 1}, 'mahalanobis') && ~strcmpi(varargin{i + 1}, 'cosine'))
                            error(strcat('Invalid parameter value (parameter #', num2str(i + 1), ').'));
                        end                        
                        obj.similarityMeasure = varargin{i + 1};
                    elseif (strcmpi(parameter, 'Echo'))
                        if (~islogical(varargin{i + 1}))
                            error(strcat('Invalid parameter value (parameter #', num2str(i + 1), ').'));
                        end
                        obj.echo = varargin{i + 1};
                    elseif (strcmpi(parameter, 'GenerateChart'))
                        if (~islogical(varargin{i + 1}))
                            error(strcat('Invalid parameter value (parameter #', num2str(i + 1), ').'));
                        end
                        obj.generateChart = varargin{i + 1}; 
                    elseif (strcmpi(parameter, 'GenerateVideo'))
                        if (~islogical(varargin{i + 1}))
                            error(strcat('Invalid parameter value (parameter #', num2str(i + 1), ').'));
                        end
                        obj.generateVideo = varargin{i + 1}; 
                    elseif (strcmpi(parameter, 'AutoMerge'))
                        if (~islogical(varargin{i + 1}))
                            error(strcat('Invalid parameter value (parameter #', num2str(i + 1), ').'));
                        end
                        obj.autoMerge = varargin{i + 1}; 
                    elseif (strcmpi(parameter, 'Delay'))
                        if (~isfloat(varargin{i + 1}))
                            error(strcat('Invalid parameter value (parameter #', num2str(i + 1), ').'));
                        end
                        obj.delay = varargin{i + 1};                     
                    end
                end
            end
        end
        
        function [threshold] = calculateThreshold(obj, s)
           if (strcmpi(obj.similarityMeasure, 'euclidean'))
               th = (obj.m ^ 2 + 1)/(2 * (s));
           elseif (strcmpi(obj.similarityMeasure, 'mahalanobis'))
               th = (obj.m ^ 2 + obj.dimension)/(2 * (s) * obj.dimension);              
           end           
           threshold = th;
        end
        
        function obj = addToChart(obj, x, y)
            if (size(x, 2) == 2)
                figure(1);                
                X = get(obj.lHandle, 'XData');
                Y = get(obj.lHandle, 'YData');
                X = [X x(1)];
                Y = [Y x(2)];
                centers = obj.getCenters(); 
                pickedColor = obj.colorPick(mod((y - 1) * 45, 256) + 1, :);
                obj.colors = [obj.colors; pickedColor];
                set(obj.lHandle, 'XData', X, 'YData', Y, 'Marker', '*', 'CData', obj.colors); 
                hold on;
                set(obj.cHandle, 'XData', centers(:, 1), 'YData', centers(:, 2), 'Marker', 'o', 'LineWidth', 5, 'MarkerEdgeColor', 'k', 'SizeData', 250);                
                drawnow;               
                hold off;
            elseif (size(x, 2) == 3)
                figure(1);                
                X = get(obj.lHandle, 'XData');
                Y = get(obj.lHandle, 'YData');
                Z = get(obj.lHandle, 'ZData');
                X = [X x(1)];
                Y = [Y x(2)];
                Z = [Z x(3)];
                centers = obj.getCenters(); 
                pickedColor = obj.colorPick(mod((y - 1) * 45, 256) + 1, :);
                obj.colors = [obj.colors; pickedColor];
                set(obj.lHandle, 'XData', X, 'YData', Y, 'ZData', Z, 'Marker', '*', 'CData', obj.colors); 
                hold on;
                set(obj.cHandle, 'XData', centers(:, 1), 'YData', centers(:, 2), 'ZData', centers(:, 3), 'Marker', 'o', 'LineWidth', 5, 'MarkerEdgeColor', 'k', 'SizeData', 250);                
                drawnow;               
                hold off;
            end
        end
        
        function [obj, y, membership] = addPoint(obj, varargin)
            nargin = size(varargin, 2);
            if (size(varargin, 2) == 0)
                error('AutoCloud needs at least one parameter.');                
            elseif (nargin > 0)
                x = varargin{1};
                if (obj.dimension == 0)
                    obj.dimension = size(x, 2);
                elseif (obj.dimension ~= size(x, 2))
                    error('The dimension of the current data sample is different from the points read so far.');
                end
                if (nargin == 2)
                    class = varargin{2};                    
                    if (ischar(class))
                        classFound = false;
                        for i = 1 : size(obj.cloudList, 2)
                            if (strcmp(obj.cloudList(i).name, class))
                                class = i;
                                classFound = true;
                                break;
                            end                                
                        end
                        if (~classFound)                            
                            error(strcat('Class "', class, '" not found.'));                        
                        end                        
                    elseif (~(floor(class) == class))
                        error(strcat('Class "', num2str(class), '" is not valid.'));
                    end  
                    while (size(obj.cloudList, 2) < class)                    
                        c = Cloud();
                        c.name = sprintf('Class %d', size(obj.cloudList, 2) + 1);
                        obj.cloudList(size(obj.cloudList, 2) + 1) = c;
                    end
                    obj.cloudList(class) = obj.cloudList(class).addPoint(x);
                    obj.k = obj.k + 1;
                    obj.intersectionList = zeros(1, size(obj.cloudList, 2));
                    obj.membershipList = zeros(1, size(obj.cloudList, 2));                    
                    obj.membershipList(class) = 1;
                    obj.intersectionList(class) = 1;
                    obj.intersectionMatrix = zeros(size(obj.cloudList, 2), size(obj.cloudList, 2));
                    y = class;
                    membership = obj.membershipList;
                    if (obj.generateChart && (obj.dimension == 2 || obj.dimension == 2))
                        obj = addToChart(x, y);
                        pause(obj.delay);
                    end
                    return;
                elseif (nargin == 1)                   
                    obj.k = obj.k + 1;
                    obj.intersectionList = zeros(1, size(obj.cloudList, 2));              
                    if(obj.k == 1)
                        c = Cloud(x);
                        c.name = 'Class 1';
                        obj.cloudList(1) = c;
                        obj.intersectionList(1) = 0;
                        obj.intersectionMatrix(1, 1) = 0;
                        obj.membershipList = 1;
                        if (obj.generateChart)
                            if (obj.dimension == 2)
                                obj.lHandle = scatter([], []);    
                                hold on;
                                obj.cHandle = scatter([], []); 
                            elseif (obj.dimension == 3)
                                obj.lHandle = scatter3([], [], []);    
                                hold on;
                                obj.cHandle = scatter3([], [], []);                             
                            else
                                warning(sprintf('Charts cannot be generated for %d-Dimensional data sets.', obj.dimension));
                            end
                        end
                    else
                        if(obj.k == 2)
                            obj.cloudList(1) = obj.cloudList(1).addPoint(x);                    
                            obj.membershipList = 1;
                        else
                             if(obj.k >= 3)
                                 createCloud = true;
                                 tauList = zeros(size(obj.cloudList,2),1);
                                 for i = 1 : size(obj.cloudList, 2)                             
                                     c1 = obj.cloudList(i);
                                     zeta = c1.calculateZeta(x, obj.similarityMeasure);                             
                                     tau = 1 - zeta;
                                     tauList(i) = tau;
                                     %%% x(k) belongs to cloud i
                                     if(zeta < Inf && zeta <= obj.calculateThreshold(c1.n))
                                         obj.cloudList(i) = obj.cloudList(i).addPoint(x);
                                         obj.intersectionList(i) = 1;
                                         createCloud = false;
                                     %%% x(k) does not belong to cloud i
                                     else
                                         obj.intersectionList(i) = 0;
                                     end
                                 end 
                                 if(createCloud == true)
                                     msg = sprintf('Creating cloud %d at instant %d.',size(obj.cloudList, 2) + 1, obj.k);
                                     if (obj.echo == true)
                                         disp(msg);
                                     end
                                     c = Cloud(x);                                     
                                     c.name = sprintf('Class %d', size(obj.cloudList, 2) + 1);                                     
                                     obj.cloudList(size(obj.cloudList, 2) + 1) = c;
                                     obj.intersectionList(size(obj.cloudList, 2)) = 1;
                                     obj.intersectionMatrix(size(obj.cloudList, 2), size(obj.cloudList, 2)) = 0;                                                          
                                 end                         
                                 obj.membershipList = tauList / sum(tauList);
                             end
                        end
                    end                    
                end
            end
            y = find(obj.membershipList == max(obj.membershipList), 1);
            if (obj.generateChart)
                obj = obj.addToChart(x, y);
                pause(obj.delay);
            end            
            if (obj.autoMerge)
                obj = mergeClouds(obj);
            end            
            membership = obj.membershipList;
        end
        
        function centers = getCenters(obj)
            centers = zeros(size(obj.cloudList, 2), obj.dimension);
            for i = 1:size(obj.cloudList, 2)
                centers(i,:) = obj.cloudList(i).mu;
            end
        end
        
        function obj = mergeClouds(obj)
            i = 1;
            while(i < size(obj.cloudList, 2))
                merge = false;                
                for j = i + 1 : size(obj.cloudList, 2)
                    if(obj.intersectionList(i) == 1 && obj.intersectionList(j) == 1)
                        obj.intersectionMatrix(i,j) = obj.intersectionMatrix(i,j) + 1;
                    end                
                    %%% recover information about clouds to be merged %%%
                    n_i = obj.cloudList(i).n;
                    n_j = obj.cloudList(j).n;
                    mu_i = obj.cloudList(i).mu;
                    mu_j = obj.cloudList(j).mu;
                    var_i = obj.cloudList(i).var;
                    var_j = obj.cloudList(j).var;
                    covmat_i = obj.cloudList(i).covmat;
                    covmat_j = obj.cloudList(j).covmat;                    
                    nint = obj.intersectionMatrix(i, j);
                    if (nint > (n_i - nint) || nint > (n_j - nint))
                        %%% merge                        
                        msg = sprintf('Merging clouds %d and %d at instant %d.', i, j, obj.k);
                        if (obj.echo == true)
                            disp(msg);
                        end
                        %%% calculate state of new cloud
                        n = n_i + n_j - nint;
                        mu = ((n_i * mu_i) + (n_j * mu_j))/(n_i + n_j);
                        var = ((n_i - 1) * var_i + (n_j - 1) * var_j)/(n_i + n_j - 2);   
                        covmat = ((n_i - 1) * covmat_i + (n_j - 1) * covmat_j)/(n_i + n_j - 2);   
                        %%% create new cloud cloud %%%
                        newCloud = Cloud();
                        newCloud = newCloud.updateCloud(mu, var, n, covmat);
                        newCloud.name = sprintf('Class %d/%d', i, j);                        
                        %%% update cloud list %%%
                        obj.cloudList = [obj.cloudList(1 : i - 1) newCloud obj.cloudList(i + 1 : j - 1) obj.cloudList(j + 1 : size(obj.cloudList, 2))];
                        
                        %%% update intersection list %%%
                        obj.intersectionList = [obj.intersectionList(1 : i - 1) 1 obj.intersectionList(i + 1 : j - 1) obj.intersectionList(j + 1 : size(obj.intersectionList, 2))];
                        %%% update intersection matrix %%%
                        A = obj.intersectionMatrix;
                        B = [A(1 : i - 1, :); zeros(1, size(A, 2)); A(i + 1 : j - 1, :); A(j + 1 : size(A, 1),:)];
                        B = [B(:, 1 : i - 1) zeros(size(B, 1), 1) B(:, i + 1 : j - 1) B(:, j + 1 : size(B, 2))];
                        C = (A(:, i) + A(:, j)).*(A(: , i).*A(:, j) ~= 0);
                        C = [C(1 : j - 1, :); C(j + 1 : size(C, 1), :)];
                        L = (A(i, :)+A(j, :)).*(A(i, :).*A(j, :) ~= 0);  
                        L = [L(:, 1 : j - 1) L(:, j + 1 : size(L, 2))];
                        B(:, i) = C;
                        B(i, :) = L;
                        B(i, i + 1 : j - 1) = A(i, i + 1 : j - 1) + A(i + 1 : j - 1, j)';                        
                        obj.intersectionMatrix = B;
                        merge = true;
                        obj.contMerge = obj.contMerge + 1;
                        break;
                    end
                end
                if(merge == true)
                    i = 1;
                else
                    i = i + 1;
                end                
            end
        end
    end
end