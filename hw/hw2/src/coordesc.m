function [b, rt, gap, tol] = coordesc(x, y, varargin)
    %COORDESC Fit a generalized coordinate descent model. The function assumes
    % that X and Y are mean centered and have standard deviation 1.
    %
    %   [B, RT, GAP, TOL] = COORDESC(X, Y) Runs coordinate descent and outputs
    %   the weights, running time, duality gap, and final tolerance.
    %
    %   Optional Parameters
    %      - alpha
    %         The L1 penalty
    %
    %      - beta
    %         The L2 penalty
    %
    %      - tol
    %         Duality gap tolerance
    %
    %      - maxiter
    %         Maximum number of iterations
    %
    %      - normalize
    %         Whether to scale the data when centering it
    
    
    % input parsing
    p = inputParser;
    p.addRequired('x', @(x) (isvector(x) || ismatrix(x)) && ~isscalar(x))
    p.addRequired('y', @(x) isvector(x))
    p.addOptional('b', zeros(size(x, 2), 1), @(x) (isempty(x) || isvector(x)) ...
        && isa(x, 'double')) 
    p.addOptional('alpha', 1, @(x) x >= 0)
    p.addOptional('beta', 0, @(x) x >= 0)
    p.addOptional('tol', 1e-3, @isposrealscalar)
    p.addOptional('maxiter', 1000, @isposintscalar)
    
    p.parse(x, y, varargin{:})
    
    r = p.Results;
    
    x = r.x;
    y = r.y;
    b = r.b;
    alpha = r.alpha;
    beta = r.beta;
    tol = r.tol;
    maxiter = r.maxiter;
    
    % get the number of samples in x (assumed to be the rows of x)
    n = size(x, 1);
            
    % we want alpha to go as high as the number of samples
    alpha = alpha * n;
    
    % setup the initial weights if not given    
    % see comments of DESCENT in descent.m for details
    tic
    [b, gap, tol] = descent(b', alpha, beta, x, y, maxiter, tol);
    rt = toc;
end
