function h = fastnewton(d)
    [X, theta, Xtheta, ~, t] = bumpgen(d);
    y = poissrnd(exp(Xtheta));
    guess = zeros(length(theta), 1);
    bb = newton(X, y, guess);
    
    hold on
    plot(t, X * bb, 'r', 'LineWidth', 2)
    plot(t, sum(X, 2), 'b--', 'LineWidth', 2)
    legend({'$F\left(\hat{\theta}\right)$', '$F\left(\theta\right)$'}, ...
        'interpreter', 'latex')
    hold off
    axis tight
    h = gcf;
end

function theta = newton(X, y, theta, varargin)
    p = inputParser;
    p.addRequired('X', @ismatrix)
    p.addRequired('y', @isvector)
    p.addRequired('theta', @isvector)
    p.addOptional('convcrit', 1e-6, @(x) isscalar(x) && isnumeric(x))
    p.addParamValue('niters', 50, @isposintscalar)
    
    p.parse(X, y, theta, varargin{:})
    r = p.Results;
    
    X = r.X;
    y = r.y;
    theta = r.theta;
    convcrit = r.convcrit;
    niters = r.niters;
    
    thetaold = 0;
    for i = 1:niters
        if sum(abs(diff(thetaold - theta))) > convcrit
            break
        end
        
        thetaold = theta;
        expXtheta = exp(X * theta);
        s = tdma(-X' * diag(expXtheta) * X, X' * (y - expXtheta));
        theta = theta(:) - s(:);
    end
end

function [X, theta, Xtheta, T, tt, ty] = bumpgen(t)
    intv = 0.01;
    k = t;
    ivals = linspace(-t, t, k);
    x = cell(1, k);
    y = cell(1, k);
    mu = zeros(1, k);
    ind = 1;
    for i = ivals
        x{ind} = i:intv:(i + 4);
        mu(ind) = mean(x{ind});
        y{ind} = normpdf(x{ind}, mu(ind), 1 / 2);
        ind = ind + 1;
    end
    x = cell2mat(x')';
    y = cell2mat(y')';
    
    tt = min(min(x)):intv:max(max(x));
    M = length(tt);
    
    i = 1;
    m = size(x, 1);
    ty = zeros(M, k);
    tx = zeros(M, k);
    offi = floor((M - m) / (k - 1));
    offmax = (offi + 1) * (k - 1);
    for offset = 1:offi:offmax
        tos = offset:(offset + (m - 1));
        tx(tos, i) = x(:, i);
        ty(tos, i) = y(:, i);
        i = i + 1;
    end
    theta = randn(k, 1);
    X = bsxfun(@times, theta', ty);
    T = tx;
    Xtheta = ty * theta;
end
