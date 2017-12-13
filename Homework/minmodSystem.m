function [result] = minmodSystem(a, b)
    n = length(a);
    result = zeros(size(a));
    for i = 1:n
        if a(i)*b(i) > 0
            if abs(a(i)) < abs(b(i))
                result(i) = a(i);
            else
                result(i) = b(i);
            end
        else
            result(i) = 0;
        end
    end
end
