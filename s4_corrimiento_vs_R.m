clear; close all; clc;
c = constantes(); p = parametros(); col = estilo(); r = rutas();

sH = escalas(c, p.ref.mL, p.ref.mN, p.ref.Z, 1);
sM = escalas(c, p.asig.mL, p.asig.mN, p.asig.Z, 1);
R = p.rp_barrido(:);
M = zeros(numel(R), 9);
for i = 1:numel(R)
    dH = corrimiento(2, 0, p.ref.Z, R(i)/sH.a, p.modelo, p.N, p.K);
    dM = corrimiento(2, 0, p.asig.Z, R(i)/sM.a, p.modelo, p.N, p.K);
    M(i, :) = [R(i), [dH.LO dH.pert dH.exacto dH.directo]*sH.Eh*1e9, ...
        [dM.LO dM.pert dM.exacto dM.directo]*sM.Eh*1e3];
    fprintf('r_p = %.2f fm listo\n', R(i));
end
cociente = (M(:, 8)*1e-3)./(M(:, 4)*1e-9);
cocienteLO = sM.coefLO(2)/sH.coefLO(2);
tabla = array2table([M cociente], 'VariableNames', {'rp_fm', 'H_LO_neV', 'H_pert_neV', ...
    'H_exacto_neV', 'H_directo_neV', 'A_LO_meV', 'A_pert_meV', 'A_exacto_meV', ...
    'A_directo_meV', 'cociente_exacto'});
writetable(tabla, fullfile(r.res, 'corrimiento_vs_R.csv'));
disp(tabla)
fprintf('Cociente de orden lider (m_r ratio)^3 Z^4 = %.6e\n', cocienteLO);

Rd = linspace(0.68, 1.02, 200);
fig = figure('Units', 'centimeters', 'Position', [2 2 26 8.5]);
subplot(1, 3, 1);
plot(Rd, sM.coefLO(2)*Rd.^2*1e3, '-', 'Color', col.gris);
hold on
plot(R, M(:, 9), 'x', 'Color', col.violeta, 'MarkerSize', 9, 'LineWidth', 1.2);
plot(R, M(:, 8), 'o', 'Color', col.naranja, 'MarkerFaceColor', col.naranja, 'MarkerSize', 5);
xlabel('r_p (fm)'); ylabel('\DeltaE_{2S} (meV)');
title(['(a) ' p.asig.nombre]);
legend({'perturbativo (ec. 8)', 'E_{ext} - E_{punt}', 'Hellmann-Feynman'}, 'Location', 'northwest');
xlim([0.68 1.02]);
subplot(1, 3, 2);
plot(Rd, sH.coefLO(2)*Rd.^2*1e9, '-', 'Color', col.gris);
hold on
plot(R, M(:, 5), 'x', 'Color', col.violeta, 'MarkerSize', 9, 'LineWidth', 1.2);
plot(R, M(:, 4), 'o', 'Color', col.azul, 'MarkerFaceColor', col.azul, 'MarkerSize', 5);
xlabel('r_p (fm)'); ylabel('\DeltaE_{2S} (neV)');
title(['(b) ' p.ref.nombre]);
legend({'perturbativo (ec. 8)', 'E_{ext} - E_{punt}', 'Hellmann-Feynman'}, 'Location', 'northwest');
xlim([0.68 1.02]);
subplot(1, 3, 3);
plot([0.68 1.02], cocienteLO*[1 1]/1e6, '-', 'Color', col.gris);
hold on
plot(R, cociente/1e6, 'o', 'Color', col.aqua, 'MarkerFaceColor', col.aqua, 'MarkerSize', 5);
xlabel('r_p (fm)'); ylabel(['\DeltaE_{2S}(' p.asig.nombre ') / \DeltaE_{2S}(' p.ref.nombre ')  (\times10^6)']);
title(['(c) Cociente entre ' char(225) 'tomos']);
legend({'(m_{r,\mu}/m_{r,e})^3', ['num' char(233) 'rico']}, 'Location', 'southwest');
xlim([0.68 1.02]);
guardar_figura(fig, 'fig3_corrimiento_vs_R');
