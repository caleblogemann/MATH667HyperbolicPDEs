function [xnew] = newton(f, df, x0, tolerance, maxIterations)
    fx = f(x0);
    xold = x0;
    mstop = 1;
    k = 0;
    while mstop && k < maxIterations
        k = k+1;
        dfx = df(xold);
        xnew = xold - fx/dfx;
        fx = f(xnew);
        if (abs(fx) < tolerance && abs(xnew - xold) < tolerance)
            mstop = 0;
        else
            xold = xnew;
        end
    end
    if (k >= maxIterations)
        disp('Newtons Method did not converge');
    end
end
