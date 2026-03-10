function x_proj = proj_box(x, a, b)
    % Projects vector x onto the box [a, b]^m
    x_proj = min(max(x, a), b);
end
