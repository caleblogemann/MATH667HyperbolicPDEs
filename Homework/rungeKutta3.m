function [sol] = rungeKutta3(L, t, w0)
    [nEqns, nCells] = size(w0);
    nTimeSteps = length(t)-1;
    deltaT = diff(t);
    k = zeros(nEqns, nCells, 3);
    sol = zeros(nEqns, nCells, nTimeSteps+1);
    sol(:,:,1) = w0;

    alpha = [1/6, 2/3, 1/6];
    lambda = [0, 0, 0; 1/2, 0, 0; -1, 2, 0];
    for n = 1:nTimeSteps
        k(:,:,1) = L(t(n), sol(:, :, n));
        k(:,:,2) = L(t(n) + sum(lambda(2,:))*deltaT(n), sol(:,:,n) + deltaT(n)*lambda(2,1)*k(:,:,1));
        k(:,:,3) = L(t(n) + sum(lambda(3,:))*deltaT(n), sol(:,:,n) + deltaT(n)*(lambda(3,1)*k(:,:,1) + lambda(3,2)*k(:,:,2)));
        sol(:,:,n+1) = sol(:,:,n) + deltaT(n)*(alpha(1)*k(:,:,1) + alpha(2)*k(:,:,2) + alpha(3)*k(:,:,3));
        disp(n/nTimeSteps);
    end

%     alpha = [1, 0, 0; 3/4, 1/4, 0; 1/3, 0, 2/3];
%     beta = [1, 0, 0; 0, 1/4, 0; 0, 0, 2/3];
%     for n = 1:nTimeSteps
%         k(:,:,1) = sol(:,:,n) + beta(1,1)*deltaT(n)*L(t(n), sol(:,:,n));
%         k(:,:,2) = alpha(2,1)*sol(:,:,n) + alpha(2,2)*k(:,:,1) + beta(2,2)*deltaT(n)*L(t(n),k(:,:,1));
%         k(:,:,3) = alpha(3,1)*sol(:,:,n) + alpha(3,3)*k(:,:,2) + beta(3,3)*deltaT(n)*L(t(n),k(:,:,2));
%         sol(:,:,n+1) = k(:,:,3);
%         disp(n/nTimeSteps);
%     end
end
