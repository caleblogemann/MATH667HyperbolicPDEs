function [u] = richtmyer(f, u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;

    a = deltaT/deltaX;
    b = 0.5 * deltaT/deltaX;

    F = zeros(nGridCells,1);
    Uhalf = zeros(nGridCells,1);

    for n = 1:nTimeSteps
        % Fluxes
        for j = 1:nGridCells
            F(j) = f(u(n, j));
        end

        % Uhalf
        % Uhalf(j) = U^(n+1/2)_(j+1/2)
        for j = 1:nGridCells-1
            Uhalf(j) = 0.5*(u(n, j)+u(n,j+1)) - b*(F(j+1) - F(j));
        end
        % zero flux boundary conditions
        Uhalf(nGridCells) = u(n, nGridCells);

        % Ustar Fluxes
        for j = 1:nGridCells-1;
            F(j) = f(Uhalf(j));
        end

        % update solution
        % zero flux boundary condition
        u(n+1, 1) = u(n, 1);
        for j = 2:nGridCells
            u(n+1,j) = u(n, j) - a*(F(j) - F(j-1));
        end
    end
end
