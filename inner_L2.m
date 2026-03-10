
% Inner product in L2 (assuming this exists)
function val = inner_L2(f, g, s)
    % Approximate L2 inner product <f,g> = integral(f*g) from 0 to 1
    val = trapz(s, f.*g);
end