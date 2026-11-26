function r = corrimiento(n, l, Z, rhoP, modelo, N, K)
g = malla(1e-4*rhoP, 20*n/Z, N);
v0 = potencial(g.rho, Z, 0, 'puntual');
dv = potencial(g.rho, Z, rhoP, modelo) - v0;
[E0, u0] = autovalor(n, l, Z, g, v0);
E1 = autovalor(n, l, Z, g, v0 + dv);
r.LO = (2/3)*Z^4*rhoP^2/n^3*(l == 0);
r.pert = trapz(g.x, u0.^2.*dv.*g.rho);
r.directo = E1 - E0;
[t, wt] = gauss_legendre(K);
r.exacto = 0;
for k = 1:K
    [~, uk] = autovalor(n, l, Z, g, v0 + t(k)*dv);
    r.exacto = r.exacto + wt(k)*trapz(g.x, uk.^2.*dv.*g.rho);
end
end
