function [E, u] = autovalor(n, l, Z, g, v)
nr = n - l - 1;
lo = -Z^2/(2*(n - 0.5)^2);
hi = -Z^2/(2*(n + 0.5)^2);
for it = 1:200
    mid = 0.5*(lo + hi);
    if mid == lo || mid == hi
        break
    end
    [~, nod] = numerov(mid, l, g, v);
    if nod > nr
        hi = mid;
    else
        lo = mid;
    end
end
E = 0.5*(lo + hi);
if nargout > 1
    w = numerov(E, l, g, v);
    u = exp(g.x/2).*w;
    u = u/sqrt(trapz(g.x, g.rho.^2.*w.^2));
end
end
