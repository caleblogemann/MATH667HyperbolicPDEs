function [L] = muscl3System(w, f, deltaX, deltaT, RFunc, LambdaFunc)
    [n, nGridCells] = size(w);
    L = zeros(n, nGridCells);
    nu = 1/deltaX;
    a = deltaX/(2*deltaT);

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

        wminus = -(1/6)*wjm1 + (5/6)*wj + (1/3)*wjp1;
        wplus = (1/3)*wj + (5/6)*wjp1 - (1/6)*wjp2;
        wtilde = wminus - wj;
        wdoubletilde = wplus + wjp1;
        wtildemod = minmod3System(wtilde, wjp1 - wj, wj - wjm1);
        wdoubletildemod = minmod3System(wdoubletilde, wjp2 - wjp1, wjp1 - wj);

        wminusmod = wj + wtildemod;
        wplusmod = wjp1 - wdoubletildemod;
        F(:,j) = a*(wminusmod - wplusmod) + 0.5*(f(wminusmod) + f(wplusmod));

        % reference solution
        %wtilde = 0.5*(wj + wjp1);
        %R = RFunc(wtilde);
        %Lambda = LambdaFunc(wtilde);
        %a = max(max(abs(Lambda)));

        %vjm1 = R\wjm1;
        %vj = R\wj;
        %vjp1 = R\wjp1;
        %vjp2 = R\wjp2;

        %vminus = -(1/6)*vjm1 + (5/6)*vj + (1/3)*vjp1;
        %vplus = (1/3)*vj + (5/6)*vjp1 - (1/6)*vjp2;
        %vtilde = vminus - vj;
        %vdoubletilde = vplus + vjp1;

        %vtildemod = minmod3System(vtilde, vjp1 - vj, vj - vjm1);
        %vdoubletildemod = minmod3System(vdoubletilde, vjp2 - vjp1, vjp1 - vj);

        %vminusmod = vj + vtildemod;
        %vplusmod = vjp1 - vdoubletildemod;

        %wminus = R*vminusmod;
        %wplus = R*vplusmod;
        %w1 = wminus - wj;
        %w2 = wplus - wj;

        %F(:,j) = f(wj + minmodSystem(w1,w2));

        %F(:,j) = a*(wminus - wplus) + 0.5*(f(wminus) + f(wplus));
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

        L(:,j) = nu*(F(:,jm1) - F(:,j));
    end
end
