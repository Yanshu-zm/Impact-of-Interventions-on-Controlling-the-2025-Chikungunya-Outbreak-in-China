function y = clip(x, a, b)
    y = min(max(x, a), b);
    %y = ceil(y);
end
