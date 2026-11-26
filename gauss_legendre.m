function [t, w] = gauss_legendre(K)
k = (1:K-1)';
b = k./sqrt(4*k.^2 - 1);
J = diag(b, 1) + diag(b, -1);
[V, D] = eig(J);
[x, idx] = sort(diag(D));
w = V(1, idx)'.^2;
t = 0.5*(x + 1);
end
