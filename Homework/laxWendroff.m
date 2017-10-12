function [u] = laxWendroff(f, df, u0, deltaT, deltaX, nTimeSteps)
    nGridCells = length(u0);
    u = zeros(nTimeSteps+1, nGridCells);
    u(1, :) = u0;
    a = 0.5 * deltaT/deltaX;
    b = 0.5 * (deltaT/deltaX)^2;

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
            % traditional laxWendroff
            Ap = df(0.5*(u(n, j) + u(n, jp1)));
            Am = df(0.5*(u(n, j) + u(n, jm1)));
            fc = f(u(n,jp1)) - f(u(n,jm1));
            fl = f(u(n,j)) - f(u(n,jm1));
            fr = f(u(n,jp1)) - f(u(n,j));
            u(n+1, j) = u(n, j) - a*fc + b*(Ap*fr - Am*fl);
        end
    end
end
