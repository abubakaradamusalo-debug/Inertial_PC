function out = lr_loss(A, b, m, x, lambda)
    Ax = A * x;
    Atran = A';
    expba = exp(-b .* Ax);
    out = Atran * (b ./ (1 + expba) - b) / m + lambda * x;
end
