clear; close all; clc;
col = estilo(); r = rutas();

T = readtable(fullfile(r.dat, 'radio_proton.csv'), 'Delimiter', ';', 'FileType', 'text');
metodos = {'e-p', 'H', 'muH', 'CODATA'};
nombres = {['dispersi' char(243) 'n e' char(8211) 'p'], ['espectroscop' char(237) 'a de H'], ['espectroscop' char(237) 'a de ' char(956) 'H'], 'ajuste CODATA'};
colores = {col.azul, col.aqua, col.naranja, col.gris};
marcas = {'s', 'd', 'o', 'p'};
desp = [-0.3, 0.3, 0, 0.6];
rotulos = {'Hand 1963', 0.35, 0, 'left'; 'Sick 2003', 0.35, 0, 'left'; 'Pohl 2010', 0.35, 0, 'left'; ...
    'Antognini 2013', 0, -0.006, 'center'; 'Bernauer 2010 (A1)', 0, 0.0115, 'center'; ...
    'Fleurbaey 2018 (1S-3S)', 0, 0.0165, 'center'; 'Xiong 2019 (PRad)', 0, -0.0185, 'center'; ...
    'Bezginov 2019 (2S-2P)', 0.35, -0.003, 'left'; 'Brandt 2022 (2S-8D)', 0.35, 0, 'left'; ...
    'Maisenbacher 2026 (2S-6P)', 0.3, 0.006, 'right'};
iA = find(strcmp(T.etiqueta, 'Antognini 2013'));
i22 = find(strcmp(T.etiqueta, 'CODATA 2022'));
i06 = find(strcmp(T.etiqueta, 'CODATA 2006'));

fig = figure('Units', 'centimeters', 'Position', [2 2 24 10]);
ax1 = axes('Position', [0.08 0.14 0.07 0.76]);
ax2 = axes('Position', [0.17 0.14 0.80 0.76]);
lim = {[1961 1965], [2001.5 2027.5]};
ejes = [ax1 ax2];
h = gobjects(1, 4);
for e = 1:2
    axes(ejes(e));
    hold on
    plot(lim{e}, T.rp_fm(i22)*[1 1], '-', 'Color', col.naranja, 'LineWidth', 0.8);
    plot(lim{e}, T.rp_fm(i06)*[1 1], '--', 'Color', col.gris, 'LineWidth', 0.8);
    for m = 1:4
        s = strcmp(T.metodo, metodos{m});
        x = T.anio(s) + desp(m);
        h(m) = errorbar(x, T.rp_fm(s), T.inc_fm(s), marcas{m}, 'Color', colores{m}, ...
            'MarkerFaceColor', colores{m}, 'MarkerSize', 5, 'LineWidth', 1.1);
    end
    for k = 1:size(rotulos, 1)
        i = find(strcmp(T.etiqueta, rotulos{k, 1}));
        m = find(strcmp(metodos, T.metodo{i}));
        xi = T.anio(i) + desp(m);
        if xi > lim{e}(1) && xi < lim{e}(2)
            text(xi + rotulos{k, 2}, T.rp_fm(i) + rotulos{k, 3}, strtrim(strtok(rotulos{k, 1}, '(')), ...
                'FontSize', 7, 'Color', col.gris, 'HorizontalAlignment', rotulos{k, 4});
        end
    end
    xlim(lim{e}); ylim([0.78 0.92]);
    box on
end
set(ax1, 'XTick', 1963);
set(ax2, 'YTickLabel', [], 'XTick', 2003:3:2027);
axes(ax1); ylabel('r_p (fm)');
axes(ax2); xlabel(['A' char(241) 'o de publicaci' char(243) 'n o de ajuste']);
text(2002, T.rp_fm(i22) - 0.004, 'CODATA 2022', 'Color', col.naranja, 'FontSize', 8);
text(2002, T.rp_fm(i06) + 0.004, 'CODATA 2006', 'Color', col.gris, 'FontSize', 8);
legend(h, nombres, 'Location', 'northeast');
title(['Historia del radio de carga del prot' char(243) 'n']);
guardar_figura(fig, 'fig6_historico_rp');

ref = T.rp_fm(iA); sref = T.inc_fm(iA);
T.tension_muH = (T.rp_fm - ref)./sqrt(T.inc_fm.^2 + sref^2);
T.tension_muH(iA) = 0;
writetable(T(:, {'anio', 'etiqueta', 'metodo', 'rp_fm', 'inc_fm', 'tension_muH'}), ...
    fullfile(r.res, 'tensiones_rp.csv'), 'Delimiter', ';');

grupos = {'H 2017-2022', strcmp(T.metodo, 'H') & T.anio >= 2017 & T.anio <= 2022; ...
    'H 2017-2026', strcmp(T.metodo, 'H') & T.anio >= 2017; ...
    'e-p 2010-2019', strcmp(T.metodo, 'e-p') & T.anio >= 2010};
S = cell(size(grupos, 1), 7);
for k = 1:size(grupos, 1)
    s = grupos{k, 2};
    w = 1./T.inc_fm(s).^2;
    mu = sum(w.*T.rp_fm(s))/sum(w);
    smu = 1/sqrt(sum(w));
    chi2 = sum(w.*(T.rp_fm(s) - mu).^2);
    nu = nnz(s) - 1;
    pval = gammainc(chi2/2, nu/2, 'upper');
    S(k, :) = {grupos{k, 1}, nnz(s), mu, smu, chi2/nu, pval, (mu - ref)/sqrt(smu^2 + sref^2)};
end
estad = cell2table(S, 'VariableNames', {'grupo', 'n', 'media_pond_fm', 'inc_fm', 'chi2_nu', 'p_valor', 'tension_vs_muH'});
writetable(estad, fullfile(r.res, 'estadistica_rp.csv'), 'Delimiter', ';');
disp(estad)

pares = {'CODATA 2006', 'Pohl 2010'; 'CODATA 2014', 'Antognini 2013'; ...
    'Bernauer 2010 (A1)', 'Xiong 2019 (PRad)'; 'Fleurbaey 2018 (1S-3S)', 'Grinin 2020 (1S-3S)'; ...
    'CODATA 2022', 'Maisenbacher 2026 (2S-6P)'};
fid = fopen(fullfile(r.res, 'tensiones_pares.txt'), 'w', 'n', 'UTF-8');
for k = 1:size(pares, 1)
    a = find(strcmp(T.etiqueta, pares{k, 1}));
    b = find(strcmp(T.etiqueta, pares{k, 2}));
    z = (T.rp_fm(a) - T.rp_fm(b))/sqrt(T.inc_fm(a)^2 + T.inc_fm(b)^2);
    for f = [1 fid]
        fprintf(f, '%-26s vs %-28s : %+.2f sigma\n', pares{k, 1}, pares{k, 2}, z);
    end
end
fclose(fid);

fid = fopen(fullfile(r.res, 'tabla_rp.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', '\begin{tabular}{rllll}', '\toprule', ...
    'A\~no & M\''etodo & $r_p$ (fm) & Tensi\''on vs.\ $\mu$H & Referencia \\', '\midrule');
mtex = containers.Map(metodos, {'Dispersi\''on $e$--$p$', 'Espectroscop\''ia H', 'Espectroscop\''ia $\mu$H', 'Ajuste CODATA'});
for i = 1:height(T)
    dec = max(0, ceil(-log10(T.inc_fm(i))) + 1);
    val = sprintf(['%.' num2str(dec) 'f(%d)'], T.rp_fm(i), round(T.inc_fm(i)*10^dec));
    if i == iA
        ten = '---';
    else
        ten = sprintf('$%+.1f\\sigma$', T.tension_muH(i));
    end
    fprintf(fid, '%d & %s & %s & %s & %s \\\\\n', T.anio(i), mtex(T.metodo{i}), val, ten, T.referencia{i});
end
fprintf(fid, '%s\n', '\bottomrule', '\end{tabular}');
fclose(fid);
