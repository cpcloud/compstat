function [X, y] = center(X, y, varargin)
    %CENTER Center and optionally scale an input matrix and vector.
    %
    %   [X, Y] = CENTER(X, Y) returns X and Y centered by their respective
    %   means, but does not scale either X or Y.
    %   
    %   Optional Arguments
    %      - normalize
    %         When this is true, X and Y will be centered by their respective
    %         means as well as scaled by their component-wise standard
    %         deviations.
    
    % variable input arguments
    p = inputParser;
    
    % X is required
    p.addRequired('X');
    
    % y is required
    p.addRequired('y');
    
    % optional normalization
    p.addOptional('normalize', false)
    
    % parse the arguments
    p.parse(X, y, varargin{:})
    
    % get the Results struct
    r = p.Results;
    
    % assign the arguments
    X = r.X;
    y = r.y;
    normalize = r.normalize;
    
    % get the mean of X
    muX = mean(X, 1);
    
    % subtract muX from each row of X
    X = bsxfun(@minus, X, muX);
    
    % if the option to normalize was given
    if normalize
        % get the standard deviation
        stdX = sqrt(sum(X .^ 2, 1));
        
        % any zeros? set them to 1
        stdX(stdX == 0) = 1;
        
        % scale X by stdX
        X = bsxfun(@rdivide, X, stdX);
    end
    
    % center y by its mean
    y = y - mean(y);
end
