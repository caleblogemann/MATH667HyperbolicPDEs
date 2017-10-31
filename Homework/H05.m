%% Problem 2 (a)
u0func = @(x) 1 + 0.5*sin(x);
du0func = @(x) 0.5*cos(x);
a = 0;
b = 2*pi;
tFinal = 1.0;
f = @(u) (u^2)/2;

E = zeros(4, 4);
iter = 0;

for N = [20, 40, 80, 160]
    iter = iter + 1;
    deltaX = (b - a)/N;
    x = linspace(a, b, N);
    u0 = u0func(x);

    deltaT = 0.5*deltaX;
    nTimeSteps = ceil(tFinal/deltaT);
    deltaT = tFinal/nTimeSteps;

    sol = godunov(f, u0, deltaT, deltaX, nTimeSteps);
    exactSol = burgersExactSolution(x, u0func, du0func, tFinal);

    E(iter, 1) = N;
    E(iter, 2) = deltaX;
    E(iter, 3) = max(abs(sol(end,:) - exactSol'));
    if(iter >= 2)
        E(iter, 4) = log(E(iter-1, 3)/E(iter, 3))/log(E(iter-1, 2)/E(iter, 2));
    end
end
disp(latexFileWriter.printMatrix(E,3));

%% Problem 2 (b)
u0func = @(x) 1 + 0.5*sin(x);
a = 0;
b = 2*pi;
tFinal = 3.0;
f = @(u) (u^2)/2;
style = ["k--", "k-"];

iter = 0;
hold on;
for N = [80,800]
    iter = iter+1;
    deltaX = (b - a)/N;
    x = linspace(a, b, N);
    u0 = u0func(x);

    deltaT = 0.5*deltaX;
    nTimeSteps = ceil(tFinal/deltaT);
    deltaT = tFinal/nTimeSteps;

    sol = godunov(f, u0, deltaT, deltaX, nTimeSteps);
    plot(x, sol(end,:), char(style(iter)), 'LineWidth', 2);
end
xlabel('x');
ylabel('u');
legend('N = 80', 'Exact Solution', 'Location', 'northwest');
hold off;
saveas(gcf, 'Figures/05_01.png', 'png');
