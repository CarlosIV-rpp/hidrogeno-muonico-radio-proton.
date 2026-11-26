function g = malla(rhoMin, rhoMax, N)
g.x = linspace(log(rhoMin), log(rhoMax), N)';
g.rho = exp(g.x);
g.h = g.x(2) - g.x(1);
end
