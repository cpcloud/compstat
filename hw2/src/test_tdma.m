function o = test_tdma(n)
    A = abs(randn(n, n));
    u = diag(A, 1);
    d = diag(A);
    l = diag(A, -1);
    x = randn(10, 1);
    A = gallery('tridiag', l, d, u);
    rhs = A * x;
    xhat = tdma([0; l], d, [u; 0], rhs);
    o = abs(x - xhat);
end
