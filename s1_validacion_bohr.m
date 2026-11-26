clear; close all; clc;
p = parametros(); col = estilo(); r = rutas();

niveles = [1 0; 2 0; 2 1; 3 0; 3 1; 3 2; 4 0; 4 3];
T = zeros(size(niveles, 1), 5);
for k = 1:size(niveles, 1)
    n = niveles(k, 1);
    l = niveles(k, 2);
    g = malla(1e-9, 30*n, p.N);
    v = potencial(g.rho, 1, 0, 'puntual');
    E = autovalor(n, l, 1, g, v);
    Eb = -1/(2*n^2);
    T(k, :) = [n, l, E, Eb, abs(E - Eb)/abs(Eb)];
end
tabla = array2table(T, 'VariableNames', {'n', 'l', 'E_num', 'E_Bohr', 'err_rel'});
writetable(tabla, fullfile(r.res, 'validacion_bohr.csv'));
disp(tabla)

etq = 'SPDFGH';
fid = fopen(fullfile(r.res, 'tabla_validacion.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', '\begin{tabular}{cccc}', '\toprule', ...
    'Estado & $E_\mathrm{num}\,[E_h]$ & $E_\mathrm{Bohr}\,[E_h]$ & Error relativo \\', '\midrule');
for k = 1:size(T, 1)
    fprintf(fid, '%d%s & $%.12f$ & $%.12f$ & %s \\\\\n', T(k, 1), etq(T(k, 2) + 1), ...
        T(k, 3), T(k, 4), tex_num(T(k, 5), 2));
end
fprintf(fid, '%s\n', '\bottomrule', '\end{tabular}');
fclose(fid);

Ns = round(logspace(log10(250), log10(8000), 11));
casos = [1 0; 2 0; 2 1];
err = zeros(numel(Ns), size(casos, 1));
for i = 1:numel(Ns)
    for k = 1:size(casos, 1)
        n = casos(k, 1);
        l = casos(k, 2);
        g = malla(1e-9, 30*n, Ns(i));
        v = potencial(g.rho, 1, 0, 'puntual');
        E = autovalor(n, l, 1, g, v);
        err(i, k) = abs(E + 1/(2*n^2))*2*n^2;
    end
end
writetable(array2table([Ns(:) err], 'VariableNames', {'N', 'err_1S', 'err_2S', 'err_2P'}), ...
    fullfile(r.res, 'convergencia.csv'));

fig = figure('Units', 'centimeters', 'Position', [2 2 14 9]);
loglog(Ns, err(:, 1), 'o-', 'Color', col.azul, 'MarkerFaceColor', col.azul, 'MarkerSize', 5);
hold on
loglog(Ns, err(:, 2), 's-', 'Color', col.naranja, 'MarkerFaceColor', col.naranja, 'MarkerSize', 5);
loglog(Ns, err(:, 3), '^-', 'Color', col.aqua, 'MarkerFaceColor', col.aqua, 'MarkerSize', 5);
loglog(Ns, err(1, 2)*(Ns/Ns(1)).^(-4), '--', 'Color', col.gris, 'LineWidth', 1);
xlabel('Puntos de malla N');
ylabel('|E_{num} - E_{Bohr}| / |E_{Bohr}|');
legend({'1S', '2S', '2P', 'N^{-4}'}, 'Location', 'southwest');
title(['Validaci' char(243) 'n con n' char(250) 'cleo puntual: Numerov en malla logar' char(237) 'tmica']);
xlim([200 1e4]);
guardar_figura(fig, 'fig1_convergencia');
