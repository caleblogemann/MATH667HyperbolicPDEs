function [u] = burgersExactSolution(x, u0, du0, t)
    nGridCells = length(x);
    u = zeros(nGridCells,1);
    g = @(x0, xj) t*u0(x0) + x0 - xj;
    dg = @(x0, xj) t*du0(x0) + 1;
    for i = 1:nGridCells
        x0 = newton(@(x0) g(x0, x(i)), @(x0) dg(x0, x(i)), x(i), 1e-5, 1000);
        u(i) = u0(x0);
    end
end
