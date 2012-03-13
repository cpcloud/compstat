function [models, x, xtheta, alphas] = test_coordesc(d, varargin)
    %TEST_COORDESC Run a test of the coordinate descent method using the
    % coordesc.m function.
    %
    %   TEST_COORDESC(d) Runs a test of the coordinate descent method by 
    %   creating a data set of d basis functions and some normally distributed 
    %   random weights. It computes the whole lasso path as it were, by 
    %   performing the coordinate descent algorithm for each value of alpha.
    %
    %   Optional Parameters
    %
    %      - nalphas
    %         The number of alphas you wish to generate
    %
    %      - epsilon
    %         Error tolerance cutoff
    %
    %      - genpdf
    %         If this is equal to true then a PDF of the resulting plot will be
    %         created in the current directory
    %
    %   Example
    %
    %      The code that created the plot in ../pdf/hw2_coordinate_descent.pdf
    %      was
    %
    %              test_coordesc(40, 1000, 'genpdf', true);
    %
    %      This says to run coordinate descent with 40 basis functions, 1000
    %      different value of alpha (montonically increasing, of course) and to
    %      generate a pdf of the resulting plot
    
    % input parsing
    ps = inputParser;
    ps.addRequired('d', @(x) x > 2 && isposintscalar(x))
    ps.addOptional('nalphas', 100, @isposintscalar)
    ps.addOptional('epsilon', 1e-3, @isposrealscalar)
    ps.addParamValue('genpdf', false, @(x) islogical(x) && isscalar(x));
    
    % get the results of parsing and put those results into variables 
    ps.parse(d, varargin{:})
    r = ps.Results;
    d = r.d;
    nalphas = r.nalphas;
    epsilon = r.epsilon;
    genpdf = r.genpdf;
   
    % generate the bump functions
    [xtheta, t, x] = bumpgen(d);
    
    % generate some data based on the bump functions
    y = poissrnd(exp(xtheta));
    
    % center the data and normalize it
    [x, y] = center(x, y, true);
    
    % weight the inputs by the outputs
    xy = x' * y;
    
    % largest alpha will be the max(abs(weighted data)) / number of samples
    alphaMax = max(abs(xy)) / size(x, 1);
    fprintf('\nalpha [min, max] == [%.4f, %.4f]\n', alphaMax * epsilon, alphaMax)
    
    % hastie et al., 2007 suggest putting the weights in a logspace
    alphas = logspace(log10(alphaMax * epsilon), log10(alphaMax), nalphas);
    
    % allocate memory for the weights matrix
    models = zeros(nalphas, d);
    
    % run coordinate descent for each alpha, putting the running time in rt
    rt = zeros(nalphas, 1);
    gap = zeros(size(rt));
    tol = zeros(size(rt));
    
    % compute the lasso path
    for i = 2:nalphas
        [models(i, :), rt(i), gap(i), tol(i)] = coordesc(x, y, ...
            models(i - 1, :), alphas(i));
    end
    
    % plot the coefficients as a function of the negative of the log a la 
    % Hastie et al. (don't remember the year...)
    interpreterOptions = {'interpreter', 'latex', 'fontsize', 15};
    xhat = x * models';
    figure
    subplot(221)
    plot(log10(alphas), models, 'linewidth', 2)
    xlabel('$\log_{10}\alpha_i$', interpreterOptions{:})
    ylabel('$\mathbf{w}$', interpreterOptions{:})
    title('Lasso', interpreterOptions{:})
    axis tight
    
    % plot the fabricated data as a function of t and as a function of the
    % middle alpha (the midpoint between extreme penalty and no penalty);
    % overlay the original data for comparison
    subplot(222)
    wxmin = xhat(:, nalphas / 2);
    hold on
    plot(t, xtheta, 'b', 'linewidth', 3)
    plot(t', wxmin, 'r', 'linewidth', 2)
    hold off
    title('``Best'''' fit', interpreterOptions{:})
    xlabel('$t$', interpreterOptions{:})
    ylabel(['$\mathbf{X}\theta$ and '...
        '$\mathbf{X}\mathbf{w}_{\mathrm{length}\left(\alpha\right)/2}$'], ...
        interpreterOptions{:})
    legend({'$\mathbf{X}\theta$', ...
        '$\mathbf{X}\mathbf{w}_{\mathrm{length}\left(\alpha\right)/2}$'}, ...
        interpreterOptions{:})
    h(2) = gca;
    
    % plot all of the models and overlay the original data
    subplot(223)
    hold on
    plot(t', xhat, 'handlevisibility', 'off');
    plot(t, xtheta, 'c', 'linewidth', 3)
    hold off
    title('Fits for each $\alpha_i$', interpreterOptions{:})
    xlabel('$t$', interpreterOptions{:})
    ylabel('$\mathbf{X}\hat\Theta$', interpreterOptions{:})
    legend({'$\mathbf{X}\theta$'}, interpreterOptions{:})

    axis tight
    h(1) = gca;
    
    linkaxes(h, 'xy')
    
    % plot the running time as a function of the sparsity of the solution
    subplot(224)
    xt = (1:length(rt))';
    [ax, h1, h2] = plotyy(xt, rt * 1000, xt, log10(alphas), @plot);
    title('Running Time and Sparsity', interpreterOptions{:})
    legend({'Running Time', '$\alpha$'}, 'interpreter', 'latex', ...
        'location', 'north')
    set(get(ax(1), 'ylabel'), 'string', 'Running Time (ms)', ...
        interpreterOptions{:})
    set(get(ax(2), 'ylabel'), 'string', '$\log_{10}\alpha_i$', ...
        interpreterOptions{:})
    set(h1, 'linewidth', 2)
    set(h2, 'linewidth', 2, 'linestyle', '--')
    
    set(ax(2), 'xtick', [])
    set(ax(2), 'xticklabel', '')
    
    % generate a PDF of the current plot, in the current directory if passed as
    % true
    if genpdf
        pdfsave('../pdf/hw2_coordinate_descent')
    end    
end
