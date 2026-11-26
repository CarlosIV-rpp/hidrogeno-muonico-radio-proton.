function [w, nodos] = numerov(E, l, g, v)
f = (l + 0.5)^2 + 2*g.rho.^2.*(v - E);
c = 1 - g.h^2*f/12;
N = numel(g.x);
w = zeros(N, 1);
w(1) = exp((l + 0.5)*(g.x(1) - g.x(2)));
w(2) = 1;
nodos = 0;
for i = 2:N-1
    w(i+1) = ((12 - 10*c(i))*w(i) - c(i-1)*w(i-1))/c(i+1);
    if w(i+1)*w(i) < 0
        nodos = nodos + 1;
    end
    if abs(w(i+1)) > 1e150
        w(1:i+1) = w(1:i+1)*1e-150;
    end
end
end
