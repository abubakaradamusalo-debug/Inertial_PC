function Fv = compute_Fv(v, t)
    % Computes A(v)(t) as defined in Example 4.3
    % v: function handle representing v(s) ∈ L²([0,1])
    % t: scalar or array of time values in [0,1]
    % Returns A(v)(t) = g(v) * F̃(v)(t)
    
    % Compute L² norm of v
    norm_v = sqrt(integral(@(s) v(s).^2, 0, 1, 'ArrayValued', true));
    
    % Compute g(v) = 1/(1 + ||v||₂²)
    g_v = 1 / (1 + norm_v^2);
    
    % Compute F̃(v)(t) = integral from 0 to t of v(s) ds
    if isscalar(t)
        F_tilde = integral(v, 0, t, 'ArrayValued', true);
    else
        F_tilde = arrayfun(@(T) integral(v, 0, T, 'ArrayValued', true), t);
    end
    
    % Compute A(v)(t) = g(v) * F̃(v)(t)
    Fv = g_v * F_tilde;
end