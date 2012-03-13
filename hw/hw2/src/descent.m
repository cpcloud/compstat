function [w, gap, tol] = descent(w, alpha, beta, X, y, maxiter, tol)
    %DESCENT Main loop of a coordinate descent solver for sparse regression.
    %
    %   [W, GAP, TOL] = descent(W, ALPHA, BETA, X, Y, MAXITER, TOL) runs the
    %   main loop of a coordinate descent solver with initial weights W, L1
    %   penalty ALPHA, L2 penalty BETA, predictors X, response variable Y. The
    %   maximum number of iterations is given by MAXITER, and the threshold for
    %   convergence is TOL. W is returned with updated weights, along with the
    %   final value of the duality gap, and final tolerance.
    
    % nfeats number of predictors (columns of X)
    nfeats = size(X, 2);
    
    % get the column norm of the data (i.e., the norm of each function)
    colnorm = sum(X .^ 2, 1);
    
    % compute the weighted residuals
    R = y - X * w;
    
    % tolerance is input tolerance * l2 norm of y
    tol = tol * norm(y) ^ 2;
    
    % init the change in tolerance per iteration
    dwtol = tol;
    
    % init the duality gap
    gap = tol + 1;
    
    % loop up to maxiter # of iterations
    for n = 1:maxiter
        % initialize the max of w and the max of the difference between 2
        % iterations of the algorithm
        wmax = 0;
        dwmax = 0;
        
        % loop over the features
        for i = 1:nfeats
            % if the column norm is zero go to the next predictor
            if ~colnorm(i)
                % this is where the algorithm takes advantage of sparsity: many
                % iterations can be skipped if the column norm of the ith 
                % feature is zero, meaning it's not really contributing that
                % much "information" to the y values
                continue
            end
            
            % store the ith weight
            wi = w(i);
            
            % if the ith weight is nonzero
            if wi
                % update the residuals with the current value + ith weight * ith
                % predictor
                R = R + wi .* X(:, i);
            end
            
            % weight the ith predictor by the updated residuals
            tmp = sum(X(:, i) .* R);
            
            % update the ith weight with the S function from Friedman et al.,
            % 2007
            w(i) = sign(tmp) * max(abs(tmp) - alpha, 0) / (colnorm(i) + beta);
            
            % if the new ith weight is nonzero
            if w(i)
                % update the residuals by subtracting the ith weight * the ith
                % predictor, i.e., penalize by w(i)
                R = R - w(i) * X(:, i);
            end
            
            % get the difference between the old weight and the current weight
            dwi = abs(w(i) - wi);
            
            % if that difference is greater than the current dwmax (recall that
            % dwmax starts at 0)
            if dwi > dwmax
                % change dwmax to dwi
                dwmax = dwi;
            end
            
            % if the ith weight is greater than the max weight (recall that wmax
            % also starts at 0)
            if abs(w(i) > wmax)
                % set wmax to the absolute value of the ith weight
                wmax = abs(w(i));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % convergence criteria %
        %%%%%%%%%%%%%%%%%%%%%%%%
        % if the max weight is nonzero OR the difference between the current
        % weights scaled by the max weight is less than the weight difference
        % threshold OR we've reached the max number of iterations
        if ~wmax || (dwmax / wmax < dwtol) || n == maxiter
            % weight the data by the residuals
            normxr = norm(X' * R - beta * w, inf);
            
            % sqrt of sum of squares of residuals
            rnorm = norm(R);
            
            % l2 norm of of w
            wnorm = norm(w);
            
            % if the weighted data is greater than the alpha penalty parameter
            if normxr > alpha
                
                % scale alpha by the weighted data
                c = alpha / normxr;
                
                % scale the residual norm with the alpha scaled by the data
                anorm = rnorm * c;
                
                % gap is the equation from the paper
                gap = 0.5 * (rnorm ^ 2 + anorm ^ 2);
            else
                % scale factor is 1
                c = 1;
                
                % gap is squared l2 norm of the residuals
                gap = rnorm ^ 2;
            end
            
            % compute the new duality gap
            gap = gap + (alpha * norm(w, 1) - c * R' * y + 0.5 * beta * ...
                (1 + c ^ 2) * wnorm ^ 2);
            
            % exit the if the dual gap is less than the specified threshold
            % before the max number of iterations is reached
            if gap < tol
                return
            end
        end
    end
    
    % exit if the gap is larger than the tolerance for the gap
    if gap > tol
        fprintf('warning: %s;\ngap = %.4f, tol = %.4f, # of iters: %n\n', ...
            'Convergence failed, try increasing the number of iterations', ...
            gap, tol, n)
    end
end
