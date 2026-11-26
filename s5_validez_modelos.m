clear; close all; clc;
c = constantes(); p = parametros(); col = estilo(); r = rutas();

rhoP = p.rhoP_barrido(:);
modelos = {'esfera', 'exponencial', 'gauss'};
Q = zeros(numel(rhoP), 4);
for i = 1:numel(rhoP)
    for m = 1:3
        d = corrimiento(2, 0, 1, rhoP(i), modelos{m}, p.N, p.K);
        Q(i, m) = d.exacto/d.LO;
        if m == 1
            Q(i, 4) = d.pert/d.LO;
        end
    end
    fprintf('r_p/a = %.2e listo\n', rhoP(i));
end
writetable(array2table([rhoP Q], 'VariableNames', {'rp_sobre_a', 'esfera_exacto', ...
    'exponencial_exacto', 'gauss_exacto', 'esfera_primer_orden'}), fullfile(r.res, 'validez_modelos.csv'));

pend = zeros(1, 3);
for m = 1:3
    sel = rhoP <= 1e-3;
    pend(m) = -(rhoP(sel)\(Q(sel, m) - 1));
end
rng(p.semilla);
Nf = 2e6;
friar = zeros(1, 3);
for m = 1:3
    d = muestra_carga(modelos{m}, Nf) - muestra_carga(modelos{m}, Nf);
    friar(m) = mean(sqrt(sum(d.^2, 2)).^3)/2;
end
fid = fopen(fullfile(r.res, 'pendientes.txt'), 'w', 'n', 'UTF-8');
for f = [1 fid]
    fprintf(f, '1 - DeltaE/DeltaE_LO ~ k (r_p/a):  k = %.3f (esfera), %.3f (exponencial), %.3f (gauss)\n', pend);
    fprintf(f, 'Friar <r^3>_(2)/(2 r_p^3) (Monte Carlo): %.3f (esfera), %.3f (exponencial), %.3f (gauss)\n', friar);
end
fclose(fid);

sH = escalas(c, p.ref.mL, p.ref.mN, p.ref.Z, p.rp);
sM = escalas(c, p.asig.mL, p.asig.mN, p.asig.Z, p.rp);
fig = figure('Units', 'centimeters', 'Position', [2 2 14 9]);
loglog(rhoP, 1 - Q(:, 1), 'o-', 'Color', col.naranja, 'MarkerFaceColor', col.naranja, 'MarkerSize', 5);
hold on
loglog(rhoP, 1 - Q(:, 2), 's-', 'Color', col.azul, 'MarkerFaceColor', col.azul, 'MarkerSize', 5);
loglog(rhoP, 1 - Q(:, 3), 'd-', 'Color', col.aqua, 'MarkerFaceColor', col.aqua, 'MarkerSize', 5);
loglog(rhoP, 1 - Q(:, 4), '--', 'Color', col.gris, 'LineWidth', 1.2);
loglog(rhoP, friar(1)*rhoP, ':', 'Color', 'k', 'LineWidth', 1);
yl = [5e-6 0.5];
plot(sH.rhoP*[1 1], yl, ':', 'Color', col.gris, 'LineWidth', 1);
plot(sM.rhoP*[1 1], yl, ':', 'Color', col.gris, 'LineWidth', 1);
text(sH.rhoP*1.2, 1.2e-5, p.ref.nombre, 'Color', col.gris);
text(sM.rhoP*1.2, 1.2e-5, p.asig.nombre, 'Color', col.gris);
set(gca, 'XScale', 'log', 'YScale', 'log');
ylim(yl);
xlabel('r_p / a');
ylabel('1 - \DeltaE_{2S} / \DeltaE_{2S}^{(ec. 8)}');
legend({'esfera uniforme', 'exponencial (dipolo)', 'gaussiana', 'esfera, solo 1er orden', ...
    ['t' char(233) 'rmino de Friar (esfera)']}, 'Location', 'northwest');
title(['Correcci' char(243) 'n m' char(225) 's all' char(225) ' de la ec. (8)']);
guardar_figura(fig, 'fig4_validez_modelos');
