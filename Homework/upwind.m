function [u] = upwind(u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;
    alpha = deltaT/deltaX;

    for n = 1:nTimeSteps
        for j = 1:nGridCells
            % periodic boundary conditions
            jm1 = j-1;
            if (j == 1)
                jm1 = nGridCells;
            end

            % update
            u(n+1, j) = u(n, j) + alpha*(u(n, jm1) - u(n, j));
        end
    end
end
