function [Xtheta, tt, ty] = bumpgen(t)
    %BUMPGEN Generate some Gaussian basis functions and some normally
    % distributed weights.
    % 
    %   XTHETA = BUMPGEN(T) Creates t functions and weights and returns the 
    %   weighted matrix of functions.
    %
    %   [XTHETA, TT] = bumpgen(t) Creates functions and weights and returns
    %   the weighted matrix and the time vector.
    %
    %   [XTHETA, TT, TY] = bumpgen(t) Creates functions and weights and returns
    %   the weighted matrix, the time vector, and the unweighted matrix.
    
    % interval
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
    Xtheta = ty * theta;
end
