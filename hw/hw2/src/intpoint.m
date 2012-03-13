function [x, rt] = intpoint(a, y, lambda, varargin)
    %INTPOINT Solve a constrained optimization problem using interior point
    % methods. This code uses a lot of the ideas in the paper included in the
    % tarball located at ../pdf/boyd-l1-regression.pdf.
    %   
    %   X = INTPOINT(A, Y, LAMBDA) returns the solution to the lasso problem.
    %
    %   [X, RT] = INTPOINT(A, Y, LAMBDA) returns the solution the lasso problem
    %   in addition to the running time of the solver.
    if ~exist(['tdma.' mexext], 'file')
        error('You need to run ''mex tdma.c'' from the command line to build the solver.')
    end
    
    % parse the input
    p = inputParser;
    
    % constraint matrix
    p.addRequired('a')
    
    % observed data
    p.addRequired('y')
    
    % penalty value
    p.addRequired('lambda')
    
    % tolerance
    p.addOptional('reltol', 1e-3)
    
    % output
    p.addOptional('x', zeros(size(a, 2), 1))
    
    % maximum newton iterations
    p.addParamValue('maxnewtoniter', 400)
    
    % max iterations for line search
    p.addParamValue('maxlinesearchiter', 200)
    
    % not sure
    p.addParamValue('mu', 2)
    
    % line search penalty values
    p.addParamValue('alpha', 0.01)
    p.addParamValue('beta', 0.5)
    
    % matrix solver to use -- defaults to the mex file produced by mex tdma.c
    p.addParamValue('solver', @tdmaimpl)
    
    p.parse(a, y, lambda, varargin{:})
    
    r = p.Results;
    a = sparse(r.a);
    y = r.y;
    lambda = r.lambda;
    reltol = r.reltol;
    x = r.x;
    
    maxnewtoniter = r.maxnewtoniter;
    maxlinesearchiter = r.maxlinesearchiter;
    mu = r.mu;
    alpha = r.alpha;
    beta = r.beta;
    solver = r.solver;
    
    % number of functions
    n = size(a, 2);
    
    % used for later internal computations
    u = ones(n, 1);
    
    % initial t
    t0 = min(max(1, 1 / lambda), 2 * n / reltol);
    t = t0;
    
    % initial f == [-1 ... -1; -1 ... -1]
    f = [x - u; -x - u];
    
    dobj = -Inf;
    s = Inf;
    
    % used later
    at = a';
    
    % initial running time
    rt = 0;
    
    % newton's method
    for i = 1:maxnewtoniter
        % compute weighted residuals
        z = a * x - y;
        
        % duality gap
        nu = 2 * z;
        maxAnu = norm(at * nu, inf);
        
        if maxAnu > lambda
            nu = nu * lambda / maxAnu;
        end
        
        % squared norm of z
        dzz = dot(z, z);
        pobj = dzz + lambda * norm(x, 1);
        dobj = max(-0.25 * dot(nu, nu) - dot(nu, y), dobj);
        gap = pobj - dobj;
        
        % stopping criterion
        if gap / dobj < reltol
            return
        end
        
        % update t
        if s >= 0.5
            t = max(min(2 * n * mu / gap, mu * t), t);
        end
        
        % newton step according to paper
        q1 = 1 ./ (u + x);
        q2 = 1 ./ (u - x);
        d1 = (q1 .^ 2 + q2 .^ 2) / t;
        d2 = (q1 .^ 2 - q2 .^ 2) / t;
        
        % gradient
        g = [at * z * 2 - (q1 - q2) / t; ...
            lambda * ones(n, 1) - (q1 + q2) / t];
        
        % make the hessian -- see paper
        dd1 = diag(d1);
        dd2 = diag(d2);
        hessian = [2 * at * a + dd1, dd2; ...
            dd2, dd1];
        
        % solve and time
        tic
        dxu = solver(hessian, -g);
        rt = rt + toc;
        
        % get the newton updates for x and u
        dx = dxu(1:n);
        du = dxu(n + 1:end);
        
        % line search
        phi = dzz + lambda * sum(u) - sum(log(-f)) / t;
        s = 1.0;
        gdx = g' * dxu;
        for j = 1:maxlinesearchiter
            newx = x + s * dx;
            newu = u + s * du;
            newf = [newx - newu; -newx - newu];
            
            if max(newf) < 0
                newz = a * newx - y;
                newphi = dot(newz, newz) + ...
                    lambda * sum(newu) - sum(log(-newf)) / t;
                
                % converged?
                if newphi - phi <= alpha * s * gdx, break, end
            end
            
            s = beta * s;
        end
        
        % if we did all of the line search iterations without convergence leave
        % newton's method loop
        if j == maxlinesearchiter, break, end
        
        x = newx;
        u = newu;
        f = newf;
    end
    
    % average the running time over the iterations
    rt = rt / i;
end

function x = tdmaimpl(A, b)
    d = b;
    a = [0; diag(A, -1)];
    b = diag(A);
    c = [diag(A, 1); 0];
    x = tdma(a, b, c, d);
end
