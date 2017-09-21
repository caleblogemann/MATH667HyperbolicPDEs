function [u] = laxWendroff(u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;
    a = 0.5 * deltaT/deltaX;
    b = 0.5 * (deltaT/deltaX)^2;

    for n = 1:nTimeSteps
        for j = 1:nGridCells
            % periodic boundary conditions
            jm1 = j-1;
            if (j == 1)
                jm1 = nGridCells;
            end
            jp1 = j+1;
            if (j == nGridCells)
                jp1 = 1;
            end

            % update
            u(n+1, j) = u(n, j) - a*(u(n,jp1) - u(n,jm1)) + b*(u(n,jp1) - 2*u(n,j) + u(n,jm1));
        end
    end
end
