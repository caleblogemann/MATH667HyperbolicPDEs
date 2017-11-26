u0func = @(x) x < 0;
Iu0func = @(x) x*(x < 0);
a = -1;
b = 2;
f = @(u) (2*u.^2)./(2*u.^2 + (1 - u).^2);
iter1 = 0;
N = 150;
tFinal = 1.0;
deltaX = (b - a)/N;
x = linspace(a+0.5*deltaX, b-0.5*deltaX, N);

u0 = zeros(N,1);
for i = 1:N
    u0(i) = (1/deltaX)*(Iu0func(x(i) + 0.5*deltaX) - Iu0func(x(i) - 0.5*deltaX));
end

deltaT = 0.5*deltaX;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;

t = 0:deltaT:tFinal;

rk3 = NumericalAnalysis.ODES.standardRK3Method;
L = @(t, u) muscl3(u, f, deltaX);
sol = rk3.solveSystem(L, t, u0);
plot(x, sol(:, end), 'k--', 'LineWidth', 2);
xlabel('x');
ylabel('u');
title(strcat('Buckley-Leverett equation at t = ', num2str(tFinal)));
saveas(gcf, 'Figures/07_01.png', 'png');
