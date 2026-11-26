clear; close all; clc;
c = constantes(); p = parametros(); col = estilo(); r = rutas();

atomos = {p.ref, p.asig};
colores = {col.azul, col.naranja};
xl = [1e-2 1e7];
yl = [1e-22 1e-1];
fig = figure('Units', 'centimeters', 'Position', [2 2 14 9]);
patch([xl(1) p.rp p.rp xl(1)], [yl(1) yl(1) yl(2) yl(2)], col.claro, 'EdgeColor', 'none');
hold on
h = gobjects(2, 1);
for k = 1:2
    A = atomos{k};
    s = escalas(c, A.mL, A.mN, A.Z, p.rp);
    g = malla(1e-9, 40/A.Z, p.N);
    v = potencial(g.rho, A.Z, 0, 'puntual');
    [~, u] = autovalor(2, 0, A.Z, g, v);
    rf = g.rho*s.a;
    P = u.^2/s.a;
    dentro = rf <= p.rp;
    Pin = trapz(g.x(dentro), u(dentro).^2.*g.rho(dentro));
    h(k) = loglog(rf, P, 'Color', colores{k}, ...
        'DisplayName', sprintf('%s:  P(r<r_p) = %.2e', A.nombre, Pin));
end
set(gca, 'XScale', 'log', 'YScale', 'log', 'Layer', 'top');
xlim(xl); ylim(yl);
text(0.7*p.rp, 1e-3, 'r < r_p', 'HorizontalAlignment', 'right', 'Color', col.gris);
xlabel('r (fm)');
ylabel('P(r) = r^2 |R_{20}(r)|^2  (fm^{-1})');
legend(h, 'Location', 'southeast');
title(sprintf('Densidad radial del estado 2S (r_p = %.5f fm)', p.rp));
guardar_figura(fig, 'fig2_densidad_radial');
