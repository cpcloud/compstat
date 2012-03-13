function rt = test_intpoint(maxd, varargin)
    %TEST_INTPOINT Run a test of the interior point solver using the intpoint.m
    % function.
    %
    %   RT = TEST_INTPOINT(maxd) Generates a data set using bumpgen(i), where
    %   i = 2 ... maxd, and then runs the interior point method on the data. RT
    %   is the running time of the solver, NOT of the whole algorithm.
    %
    %   Optional parameters
    %
    %      - lambda
    %           A scalar or vector of values of the penalty that you want to
    %           give to the solver.
    %
    %      - genpdf
    %           If this is true, then the script will generate a PDF of the 
    %           resulting plot in the current directory.
    %   
    %   Example
    %
    %      The code to generate the plots in ../pdf/hw2_interior_point.pdf is
    %      
    %                    test_intpoint(40, 'genpdf', true)
    %
    %      This runs the interior point method for sparse regression 3 times,
    %      because there are by default 3 values in the lambda vector, with 40
    %      basis functions. It generates a PDF of the resulting plot.
    
    p = inputParser;
    p.addRequired('maxd')
    p.addOptional('lambda', [10, 50, 100])
    p.addParamValue('genpdf', false, @(x) islogical(x) && isscalar(x) || ...
        (isinteger(x) && (x == 1 || x == 0)));
    p.parse(maxd, varargin{:})
    r = p.Results;
    maxd = r.maxd;
    lambda = r.lambda;
    genpdf = r.genpdf;
    
    rt = zeros(maxd, length(lambda));
    
    % for each d
    for d = 2:maxd
        % generate d functions and some weights
        [xtheta, t, x] = bumpgen(d);
        
        % for each lambda
        for j = 1:length(lambda)
            % run the interior point solver
            [xhat{j}, rt(d - 1, j)] = intpoint(x, xtheta, lambda(j));
        end
    end
    
    % get the average running time over lambdas
    rt = mean(rt, 2);
    
    % plotting
    figure
    for i = 1:length(lambda)
        subplot(2, 3, i)
        hold on
        
        % show the 'true' weighted data vs each fit 
        plot(t', xtheta, 'b', 'linewidth', 5)
        plot(t', x * xhat{i}, 'r', 'linewidth', 3)
        
        hold off
        
        legend({'$\mathbf{X}\theta$', '$\mathbf{X}\hat{\theta}$'}, ...
            'interpreter', 'latex', 'fontsize', 10)
        xlabel('$t$', 'interpreter', 'latex', 'fontsize', 10)
        ylabel('$\mathbf{X}\theta$', 'interpreter', 'latex', 'fontsize', 10)
        title(['$\lambda=' num2str(lambda(i)) ',\,\,d=' num2str(maxd) '$'], ...
            'interpreter', 'latex', 'fontsize', 10)
        axis tight
    end
    
    % plotting the running time of the solver
    subplot(2, 3, 4:6)
    plot(rt * 1e3, 'r', 'linewidth', 2)
    xlabel('$d$', 'interpreter', 'latex', 'fontsize', 15)
    ylabel('Running Time (ms)', 'interpreter', 'latex', 'fontsize', 15)
    title('C Code Solver', 'interpreter', 'latex', 'fontsize', 15)
    
    % if the parameter genpdf is passed as true then make a plot
    if genpdf
        pdfsave('../pdf/hw2_interior_point')
    end
end
