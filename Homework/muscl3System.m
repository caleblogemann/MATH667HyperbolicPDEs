function [L] = muscl3System(w, f, deltaX, RFunc, LambdaFunc)
    [n, nGridCells] = size(w);
    L = zeros(n, nGridCells);
    nu = 1/deltaX;

    boundaryConditions = 'zeroFlux';

    % flux array, F(:, i) is flux at i + 1/2 interface
    F = zeros(n, nGridCells);

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
        wjm1 = w(:,jm1);
        wj = w(:,j);
        wjp1 = w(:,jp1);
        wjp2 = w(:,jp2);

        % reference solution
        wtilde = 0.5*(wj + wjp1);
        R = RFunc(wtilde);
        Lambda = LambdaFunc(wtilde);

        vjm1 = R\wjm1;
        vj = R\wj;
        vjp1 = R\wjp1;
        vjp2 = R\wjp2;

        vminus = -(1/6)*vjm1 + (5/6)*vj + (1/3)*vjp1;
        vplus = (1/3)*vj + (5/6)*vjp1 - (1/6)*vjp2;
        vtilde = vminus - vj;
        vdoubletilde = vplus + vjp1;

        vtildemod = minmod3(vtilde, vjp1 - vj, vj - vjm1);
        vdoubletildemod = minmod3(vdoubletilde, vjp2 - vjp1, vjp1 - vj);

        vminusmod = vj + vtildemod;
        vplusmod = vjp1 - vdoubletildemod;

        wminus = R*vminusmod;
        wplus = R*vplusmod;

        v1 = vminusmod - vj;
        v2 = vplusmod - vj;

        F(j) = f(vj + minmod(v1, v2));
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
