u0func = @(x) 1.0*(x <= 1.0);
deltaX = 0.001;
a = 0.0;
b = 2.0;
nGridCells = (b - a)/deltaX;

deltaT = deltaX/2.0;
x = linspace(a, b, nGridCells);
u0 = u0func(x);

tFinal = 5;
nTimeSteps = tFinal/deltaT;

uExactFunc = @(x, t) 1.0*(mod(x - t, 2) <= 1.0);
upwindSol = upwind(u0, deltaT, deltaX, nTimeSteps);
laxFriedrichsSol = laxFriedrichs(u0, deltaT, deltaX, nTimeSteps);
laxWendroffSol = laxWendroff(u0, deltaT, deltaX, nTimeSteps);

t2 = 2/deltaT;
t5 = nTimeSteps;

uExactSol2 = arrayfun(@(xj) uExactFunc(xj, 2), x);
plot(x, upwindSol(t2+1,:), 'k--', x, laxFriedrichsSol(t2+1,:), 'k+', ...
    x,  laxWendroffSol(t2+1,:), 'ko', x, uExactSol2, 'k-');
xlabel('x');
ylabel('u');
title('T = 2');
legend('Upwind', 'Lax Friedrichs', 'Lax Wendroff', 'Exact');
saveas(gcf, 'Figures/02_01.png', 'png');

uExactSol5 = arrayfun(@(xj) uExactFunc(xj, 5), x);
plot(x, upwindSol(t5+1,:), 'k--', x, laxFriedrichsSol(t5+1,:), 'k+', ...
    x, laxWendroffSol(t5+1,:), 'ko', x, uExactSol5, 'k-');
xlabel('x');
ylabel('u');
title('T = 5');
legend('Upwind', 'Lax Friedrichs', 'Lax Wendroff', 'Exact');
saveas(gcf, 'Figures/02_02.png', 'png');

