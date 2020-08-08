% proximal operator of r(x) = y*|x|_1

function Z = shrink_matrix(X, y)
n = length(X);
Z = zeros(n, n);
for i = 1:n
    for j = 1:n
        if X(i,j) > y
            Z(i,j) = X(i,j) - y;
        end
        if X(i,j) < -y
            Z(i,j) = X(i,j) + y;
        end
    end
end 