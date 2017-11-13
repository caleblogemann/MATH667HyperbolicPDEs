function [u] = centralFV2(f, u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;
    nu = deltaT/deltaX;

    boundaryConditions = 'periodic';

    % flux array, F(i) is flux at i - 1/2 interface
    F = zeros(nGridCells,1);

    for n = 1:nTimeSteps
        % compute fluxes at boundaries
        for j = 1:nGridCells
            % boundary conditions
            jm1 = j-1;
            if (j == 1)
                if (strcmp(boundaryConditions,'periodic'))
                    jm1 = nGridCells;
                elseif (strcmp(boundaryConditions,'zeroFlux'))
                    jm1 = 1;
                end
            end
            F(j) = f(0.5*(u(n,jm1) + u(n,j)));
        end

        % update solution
        for j = 1:nGridCells
            % boundary conditions
            jp1 = j+1;
            if (j == nGridCells)
                if (strcmp(boundaryConditions,'periodic'))
                    jp1 = 1;
                elseif (strcmp(boundaryConditions,'zeroFlux'))
                    jp1 = nGridCells;
                end
            end

            u(n+1, j) = u(n, j) + nu*(F(j) - F(jp1));
        end
    end
end
