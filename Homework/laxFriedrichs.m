function [u] = laxFriedrichs(f, u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;
    a = 0.5 * deltaT/deltaX;

    for n = 1:nTimeSteps
        for j = 1:nGridCells
            % zero flux boundary conditions
            jm1 = j-1;
            if (j == 1)
                jm1 = 1;
            end
            jp1 = j+1;
            if (j == nGridCells)
                jp1 = nGridCells;
            end

            % update
            % u(n+1, j) = u(n, j) - a*(u(n,jp1) - u(n,jm1)) + 0.5*(u(n,jp1) - 2*u(n,j) + u(n,jm1));
            u(n+1, j) = 0.5*(u(n,jm1) + u(n,jp1)) - a*(f(u(n,jp1)) - f(u(n,jm1)));
        end
    end
end
