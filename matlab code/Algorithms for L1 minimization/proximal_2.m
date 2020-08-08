% proximal operator of r(x) = y*|x|_2

function z = proximal_2(x, y)
n = length(x);
z = zeros(n, 1);
m = norm(x, 2);
if m > y
    z = ((m - y) / m) * x;
end