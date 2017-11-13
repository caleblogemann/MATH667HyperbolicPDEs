function [result] = minmod(a, b)
    if a*b > 0
        if abs(a) < abs(b)
            result = a;
        else
            result = b;
        end
    else
        result = 0;
    end
end
