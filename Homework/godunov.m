function [u] = godunov(f, u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;
    nu = deltaT/deltaX;

    % flux array, F(i) is flux at i - 1/2 interface
    F = zeros(nGridCells,1);

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
            if (u(n, jm1) <= u(n, j))
                % min[ul < u < ur]{f(u)}
                % specific to burger's equation
                if (u(n, jm1)*u(n, j) < 0)
                    F(j) = 0;
                else
                    F(j) = min(ful, fur);
                end
            else
                % max[ur < u < ul]{f(u)}
                F(j) = max(ful, fur);
            end
        end

        % update solution
        for j = 1:nGridCells
            % zero flux boundary conditions
            jp1 = j+1;
            if (j == nGridCells)
                jp1 = nGridCells;
            end

            u(n+1, j) = u(n, j) + nu*(F(j) - F(jp1));
        end
    end
end
