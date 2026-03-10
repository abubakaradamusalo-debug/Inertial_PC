function Fx = F(x)
    m = length(x);
    
    % Enforce boundary conditions
    x_ext = [0; x(:); 0];  % x_0 = 0, x_{m+1} = 0

    % F1(x) computation
    F1 = zeros(m, 1);
    for i = 1:m
        xi_m1 = x_ext(i);     % x_{i-1}
        xi    = x_ext(i+1);   % x_i
        xi_p1 = x_ext(i+2);   % x_{i+1}
        F1(i) = xi_m1^2 + xi^2 + xi_m1*xi + xi*xi_p1;
    end

    % F2(x) = D*x + d
    D = zeros(m, m);
    for i = 1:m
        for j = 1:m
            if i == j
                D(i,j) = 4;
            elseif i - j == 1
                D(i,j) = 1;
            elseif i - j == -1
                D(i,j) = -2;
            else
                D(i,j) = 0;
            end
        end
    end

    d = -ones(m, 1);
    F2 = D * x + d;

    % Final operator
    Fx = F1 + F2;
end
