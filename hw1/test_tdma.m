function m = test_tdma(n, varargin)
    p = inputParser;
    p.addRequired('n', @(x) isposintscalar(x) && x > 1)
    p.addOptional('k', 100, @(x) isposintscalar(x) && x > 1)
    p.parse(n, varargin{:})
    r = p.Results;
    n = r.n;
    k = r.k;
    
    h = zeros(n, k);
    for i = 1:n
        h(i, :) = test_tdma_impl(k);
    end
    m = median(h, 1);
end

function t = test_tdma_impl(d, varargin)
    p = inputParser;
    p.addRequired('d', @(x) isposintscalar(x) && x > 1)
    p.addOptional('solver', @tdma, @(x) isa(x, 'function_handle'));
    p.addOptional('randmax', 10, @(x) x > 2);
    p.parse(d, varargin{:})
    r = p.Results;
    d = r.d;
    randmax = r.randmax;
    solver = r.solver;
    t = zeros(d, 1);
    parfor i = 2:d
        l = randn + zeros(i - 1, 1);
        u = l;
        m = abs(gallery('tridiag', i, l, randn(i, 1), u));
        b = randi(randmax, [size(m, 1), 1]);
        t(i) = timeit(@() solver(m, b)); %#ok<PFBNS>
    end
end
