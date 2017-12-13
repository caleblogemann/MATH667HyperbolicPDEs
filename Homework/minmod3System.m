function [result] = minmod3System(a, b, c)
    n = length(a);
    result = zeros(size(a));
    for i = 1:n
        temp = [a(i), b(i), c(i)];
        if a(i)*b(i) > 0 && b(i)*c(i) > 0
            [~,idx] = min(abs(temp));
            result(i) = temp(idx);
        else
            result(i) = 0;
        end
    end
end
