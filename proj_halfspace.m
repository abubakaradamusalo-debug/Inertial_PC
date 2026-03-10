function pj = proj_halfspace(x, n, y)
    % Function to project a point x onto the halfspace defined by
    % <n, x - y> <= 0, where n is the normal vector and y is a point on the
    % boundary of the halfspace.
    %
    % Inputs:
    %   x - point to be projected (vector)
    %   n - normal vector defining the halfspace (vector)
    %   y - point on the boundary of the halfspace (vector)
    %
    % Output:
    %   pj - projection of point x onto the halfspace

    % Check if x is already in the halfspace
    if dot(n, x - y) <= 0
        pj = x;  % x is already in the halfspace
    else
        % Project onto the halfspace
        pj = x - dot(n, x - y) / norm(n)^2 * n;
    end
end
