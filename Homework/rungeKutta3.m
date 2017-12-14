function [sol] = rungeKutta3(L, t, w0)
    alpha = [1/6, 2/3, 1/6];
    lambda = [0, 0, 0; 1/2, 0, 0; -1, 2, 0];

    [nEqns, nCells] = size(w0);
    nTimeSteps = length(t);
    deltaT = diff(t);
    k = zeros(nEqns, nCells, 3);

    for n = 1:nTimeSteps
        k(:,:,1) = L(t(n), sol(:, :, n))
        k(:,:,2) = L(t(n) + sum(lambda(2,:))*deltaT(n), sol(:,:,n) + deltaT(n)*lambda(2,1)*k(:,:,1));
        k(:,:,2) = L(t(n) + sum(lambda(3,:))*deltaT(n), sol(:,:,n) + deltaT(n)*(lambda(3,1)*k(:,:,1) + lambda(3,2)*k(:,:,2)));
        sol(:,:,n+1) = sol(:,:,n) + deltaT(n)*(alpha(1)*k(:,:,1) + alpha(2)*k(:,:,2) + alpha(3)*k(:,:,3));
    end
end
