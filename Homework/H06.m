%% Problem 1 (a)
u0func = @(x) 1 + 0.5*sin(x);
Iu0func = @(x) x - 0.5*cos(x);
du0func = @(x) 0.5*cos(x);
a = 0;
b = 2*pi;
tFinal = 1.0;
f = @(u) (u^2)/2;

for method = ["centralFV2", "upwindFV2", "muscl2"]
    E = zeros(4, 4);
    iter = 0;

    for N = [20, 40, 80, 160]
        iter = iter + 1;
        deltaX = (b - a)/N;
        x = linspace(a+0.5*deltaX, b-0.5*deltaX, N);
        u0 = zeros(N,1);
        for i = 1:N
            u0(i) = (1/deltaX)*(Iu0func(x(i) + 0.5*deltaX) - Iu0func(x(i) - 0.5*deltaX));
        end

        deltaT = 0.25*deltaX;
        nTimeSteps = ceil(tFinal/deltaT);
        deltaT = tFinal/nTimeSteps;

        sol = feval(char(method), f, u0, deltaT, deltaX, nTimeSteps);
        exactSol = burgersExactSolution(x, u0func, du0func, tFinal);
        plot(x, sol(end,:), x, exactSol);
        pause();
        

        E(iter, 1) = N;
        E(iter, 2) = deltaX;
        E(iter, 3) = max(abs(sol(end,:) - exactSol'));
        if(iter >= 2)
            E(iter, 4) = log(E(iter-1, 3)/E(iter, 3))/log(E(iter-1, 2)/E(iter, 2));
        end
    end
    disp(latexFileWriter.printMatrix(E,3));
end

%% Problem 1 (b)
u0func = @(x) 1 + 0.5*sin(x);
Iu0func = @(x) x - 0.5*cos(x);
a = 0;
b = 2*pi;
tFinal = 3.0;
f = @(u) (u^2)/2;
style = ["k--", "k-"];

% exact solution
N = 800;
deltaX = (b - a)/N;
xExact = linspace(a, b, N);
u0 = u0func(xExact);

deltaT = 0.5*deltaX;
nTimeSteps = ceil(tFinal/deltaT);
deltaT = tFinal/nTimeSteps;
exactSol = godunov(f, u0, deltaT, deltaX, nTimeSteps);

iter1 = 0;
for method = ["centralFV2", "upwindFV2", "muscl2"]
    iter1 = iter1+1;
    iter = 0;
    N = 80;
    iter = iter+1;
    deltaX = (b - a)/N;
    x = linspace(a+0.5*deltaX, b-0.5*deltaX, N);
    u0 = zeros(N,1);
    for i = 1:N
        u0(i) = (1/deltaX)*(Iu0func(x(i) + 0.5*deltaX) - Iu0func(x(i) - 0.5*deltaX));
    end

    deltaT = 0.1*deltaX;
    nTimeSteps = ceil(tFinal/deltaT);
    deltaT = tFinal/nTimeSteps;

    sol = feval(char(method), f, u0, deltaT, deltaX, nTimeSteps);
    plot(x, sol(end,:), 'k--', xExact, exactSol(end,:), 'k-', 'LineWidth', 2);
    xlabel('x');
    ylabel('u');
    title(char(method));
    legend('N = 80', 'Exact Solution', 'Location', 'northwest');
    saveas(gcf, strcat('Figures/06_0',num2str(iter1),'.png'), 'png');
end

%% Problem 2
%u0func = @(x) exp(-x.^2);
%Iu0func = @(x) sqrt(pi)/2*erf(x);
u0func = @(x) (x.^2)*(-1 < x && x < 1);
Iu0func = @(x) (1/3*x.^3)*(-1 < x && x < 1);
a = -pi;
b = pi;
f = @(u) u;
exactSol = @(x, t) u0func(x - t);
style = ["k--", "k-"];
iter1 = 0;
for tFinal = [1.0, 2.0]
    for method = ["upwindFV2", "muscl2"]
        iter1 = iter1+1;
        figure;
        hold on;
        N = 80;
        iter = iter+1;
        deltaX = (b - a)/N;
        x = linspace(a+0.5*deltaX, b-0.5*deltaX, N);
        u0 = zeros(N,1);
        for i = 1:N
            u0(i) = (1/deltaX)*(Iu0func(x(i) + 0.5*deltaX) - Iu0func(x(i) - 0.5*deltaX));
        end

        deltaT = 0.01*deltaX;
        nTimeSteps = ceil(tFinal/deltaT);
        deltaT = tFinal/nTimeSteps;

        sol = feval(char(method),f, u0, deltaT, deltaX, nTimeSteps);
        plot(x, sol(end,:), 'k--', 'LineWidth', 2);
        exactSolution = zeros(N,1);
        for i = 1:N
            exactSolution(i) = exactSol(x(i), tFinal);
        end
        plot(x, exactSolution, 'k-', 'LineWidth', 2);
        xlabel('x');
        ylabel('u');
        title(strcat(char(method),' at t = ', num2str(tFinal)));
        legend('N = 80', 'Exact Solution', 'Location', 'northwest');
        hold off;
        saveas(gcf, strcat('Figures/06_0',num2str(iter1+3),'.png'), 'png');
    end
end
