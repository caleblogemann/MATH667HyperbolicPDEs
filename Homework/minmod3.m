function [result] = minmod3(a, b, c)
    temp = [a, b, c];
    if a*b > 0 && b*c > 0
        [~,idx] = min(abs(temp));
        result = temp(idx);
    else
        result = 0;
    end
end
