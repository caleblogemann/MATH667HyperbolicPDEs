function [u] = roe(f, u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;
    nu = deltaT/deltaX;

    % flux array, F(i) is flux at i - 1/2 interface
    % periodic boundary conditions F(nGridCells + 1) = F(1), so F(nGridCells + 1) not necessary
    F = zeros(nGridCells);

    for n = 1:nTimeSteps
        % compute fluxes at boundaries
        for j = 1:nGridCells
            % zero flux boundary conditions
            jm1 = j-1;
            if (j == 1)
                jm1 = 1;
            end

            ful = f(u(n,jm1));
            fur = f(u(n, j));
            if ((ful - fur)/(u(n,jm1) - u(n,j)) >= 0)
                F(j) = ful;
            else
                F(j) = fur;
            end
        end

        % update solution
        for j = 1:nGridCells
            % periodic boundary conditions
            jp1 = j+1;
            if (j == nGridCells)
                jp1 = nGridCells;
            end

            u(n+1, j) = u(n, j) + nu*(F(j) - F(jp1));
        end
    end
end
