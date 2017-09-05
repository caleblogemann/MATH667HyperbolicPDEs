u0func = @(x) 1.0*(x <= 1.0);
deltaX = 0.01;
a = 0.0;
b = 2.0;
nGridCells = (b - a)/deltaX;

deltaT = deltaX/2.0;
x = linspace(a, b, nGridCells);
u0 = u0func(x);

tFinal = 50;
nTimeSteps = tFinal/deltaT;
u = upwind(u0, deltaT, deltaX, nTimeSteps);

t2 = 2/deltaT;
t10 = 10/deltaT;
t50 = nTimeSteps;
plot(x, u(t2+1,:), x, u(t10+1,:), x, u(t50+1,:));
legend('T = 2', 'T = 10', 'T = 50');
xlabel('x');
ylabel('u');
title('Forward Euler Upwind');
saveas(gcf, 'Figures/01_02.png', 'png');
