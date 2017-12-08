function [L] = muscl3System(u, f, deltaX)
    nGridCells = length(u);
    L = zeros(nGridCells, 1);
    nu = 1/deltaX;

    boundaryConditions = 'zeroFlux';

    % flux array, F(i) is flux at i + 1/2 interface
    F = zeros(nGridCells,1);

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

        jp1 = j+1;
        jp2 = j+2;
        if (j == nGridCells)
            if (strcmp(boundaryConditions,'periodic'))
                jp1 = 1;
                jp2 = 2;
            elseif (strcmp(boundaryConditions,'zeroFlux'))
                jp1 = nGridCells;
                jp2 = nGridCells;
            end
        elseif (j == nGridCells - 1)
            if (strcmp(boundaryConditions,'periodic'))
                jp2 = 1;
            elseif (strcmp(boundaryConditions,'zeroFlux'))
                jp2 = nGridCells;
            end
        end

        uminus = -(1/6)*u(jm1) + (5/6)*u(j) + (1/3)*u(jp1);
        uplus = (1/3)*u(j) + (5/6)*u(jp1) - (1/6)*u(jp2);
        utilde = uminus - u(j);
        udoubletilde = uplus + u(jp1);

        utildemod = minmod3(utilde, u(jp1) - u(j), u(j) - u(jm1));
        udoubletildemod = minmod3(udoubletilde, u(jp2) - u(jp1), u(jp1) - u(j));

        uminusmod = u(j) + utildemod;
        uplusmod = u(jp1) - udoubletildemod;

        u1 = uminusmod - u(j);
        u2 = uplusmod - u(j);
        F(j) = f(u(j) + minmod(u1, u2));
    end

    % update solution
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

        L(j) = nu*(F(jm1) - F(j));
    end
end
