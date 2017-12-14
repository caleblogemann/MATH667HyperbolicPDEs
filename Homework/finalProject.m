g = 1.4;
energyFuncPrimitive = @(rho, u, P) P/(g - 1) + 0.5*rho*u^2;
energyFunc = @(rho, m, P) P/(g - 1) + 0.5*m^2/rho;
pressureFuncPrimitive = @(rho, u, E) (g - 1)*(E - 0.5*rho*u^2);
pressureFunc = @(rho, m, E) (g - 1)*(E - 0.5*m^2/rho);

rhoL = 1;
uL = 0;
mL = rhoL*uL;
PL = 1;
EL = energyFuncPrimitive(rhoL, uL, PL);
rhoR = 0.125;
uR = 0;
mR = rhoR*uR;
PR = 0.1;
ER = energyFuncPrimitive(rhoR, uR, PR);

w0func = @(x) [rhoL, mL, EL]'*(x <= 0) + [rhoR, mR, ER]'*(x > 0);
a = -5;
b = 5;
fConservative = @(rho, m, E) [m; m^2/rho + pressureFunc(rho, m, E); m/rho*(E + pressureFunc(rho, m, E))];
fPrimitive = @(rho, u, P) [rho*u; rho*u^2 + P; u*(energyFuncPrimitive(rho, u, P) + P)];
f = @(w) fConservative(w(1), w(2), w(3));
N = 1000;
tFinal = 2.0;
deltaX = (b - a)/N;
x = linspace(a+0.5*deltaX, b-0.5*deltaX, N);
w0 = w0func(x);

RPrimitive = @(rho, u, P) [1, 1, 1; u - sqrt(g*P/rho), u, u + sqrt(g*P/rho); g/(g - 1)*P/rho + u^2/2 - u*sqrt(g*P/rho), 0.5*u^2, g/(g - 1)*P/rho + u^2/2 - u*sqrt(g*P/rho)];
RConservative = @(rho, m, E) RPrimitive(rho, m/rho, pressureFunc(rho, m, E));
RFunc = @(w) RConservative(w(1), w(2), w(3));
multByR = @(w) RFunc(w)*w;
%multByRInverse = @(w) RFunc(w)\w;

LambdaPrimitive = @(rho, u, P) [u - sqrt(g*P/rho), 0, 0; 0, u, 0; 0, 0, u + sqrt(g*P/rho)];
LambdaConservative = @(rho, m, E) LambdaPrimitive(rho, m/rho, pressureFunc(rho, m, E));
LambdaFunc = @(w) LambdaConservative(w(1), w(2), w(3));
multByLambda = @(w) LambdaFunc(w)*w;

cfl = 0.1;
deltaT = cfl*deltaX;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;

t = 0:deltaT:tFinal;

%v0 = multByRInverse(w0);

L = @(t, u) muscl3System(u, f, deltaX, deltaT, RFunc);
%rk3 = NumericalAnalysis.ODES.standardRK3Method;
%sol = rk3.solveSystem(L, t, w0);
sol = rungeKutta3(L, t, w0);
rho = sol(1,:,end);
u = sol(2,:,end)./sol(1,:,end);
p = (g - 1)*(sol(3,:,end) - 0.5*rho.*u.^2);
subplot(3, 1, 1);
plot(x, rho, 'LineWidth', 2);
xlabel('x');
ylabel('rho');
title('Density');
subplot(3, 1, 2);
plot(x, u, 'LineWidth', 2);
xlabel('x');
ylabel('u');
title('Velocity');
subplot(3, 1, 3);
plot(x, p, 'LineWidth', 2);
xlabel('x');
ylabel('P');
title('Pressure');
saveas(gcf, 'Figures/finalProject.png', 'png');
figure;
plot(x, sol(:, :,end), 'k--', 'LineWidth', 2);
xlabel('x');
ylabel('w');
title(strcat('Euler equations at t = ', num2str(tFinal)));
saveas(gcf, 'Figures/finalProject_2.png', 'png');
