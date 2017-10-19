function [u] = localLaxFriedrichs(f, df, u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;
    a = deltaT/deltaX;

    F = zeros(nGridCells);

    for n = 1:nTimeSteps
        % F(j) = f_{j-1/2}
        for j = 1:nGridCells
            % zero flux boundary conditions
            jm1 = j-1;
            if (j == 1)
                jm1 = 1;
            end
            alpha = max(abs(df(u(n, jm1))), abs(df(u(n, j))));
            F(j) = 0.5*(f(u(n,jm1)) + f(u(n, j)) - alpha*(u(n, j) - u(n, jm1)));
        end

        for j = 1:nGridCells
            % zero flux boundary conditions
            jp1 = j+1;
            if (j == nGridCells)
                jp1 = nGridCells;
            end

            % update
            % u(n+1, j) = u(n, j) - a*(u(n,jp1) - u(n,jm1)) + 0.5*(u(n,jp1) - 2*u(n,j) + u(n,jm1));
            u(n+1, j) = u(n, j) - a*(F(jp1) - F(j));
        end
    end
end
