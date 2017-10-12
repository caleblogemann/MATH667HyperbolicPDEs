%% 2 (a)
f = @(u) (0.5*u^2);
df = @(u) u;
%u0func = @(x) 1.0*(x <= 0.0) - 0.5*(x > 0.0);
u0func = @(x) -1.0*(x < 0.0) + 0.5*(x > 0.0);
deltaX = 0.01;
a = -5.0;
b = 5.0;
nGridCells = (b - a)/deltaX;

deltaT = deltaX/2.0;
%x = linspace(a+0.5*deltaX, b-0.5*deltaX, nGridCells-1);
x = linspace(a, b, nGridCells);
u0 = u0func(x);

tFinal = 2;
nTimeSteps = tFinal/deltaT;

uExactFunc = @(x, t) -1.0*(x <= -t) + (x/t)*(x > -t && x < 0.5*t) + 0.5*(x >= 0.5*t);
%uExactFunc = @(x, t) 1.0*(x <= 0.25*t) - 0.5*(x >= 0.25*t);
roeSol = roe(f, u0, deltaT, deltaX, nTimeSteps);
laxFriedrichsSol = laxFriedrichs(f, u0, deltaT, deltaX, nTimeSteps);
laxWendroffSol = richtmyer(f, u0, deltaT, deltaX, nTimeSteps);

uExactSol = arrayfun(@(xj) uExactFunc(xj, tFinal), x);
plot(x, roeSol(nTimeSteps+1,:), 'k--', x, laxFriedrichsSol(nTimeSteps+1,:), 'k-.', x, laxWendroffSol(nTimeSteps+1,:), 'k:', x, uExactSol, 'k-','LineWidth',2);
xlabel('x');
ylabel('u');
%xlim([0.4, 0.6]);
xlim([-2.5, 1.5]);
ylim([-1.1, 0.6]);
title('T = 2');
legend('Roe', 'Lax-Friedrichs', 'Richtmyer', 'Exact', 'Location', 'southeast');
saveas(gcf, 'Figures/03_02.png', 'png');
