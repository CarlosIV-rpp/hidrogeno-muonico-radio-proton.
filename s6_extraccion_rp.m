clear; close all; clc;
c = constantes(); p = parametros(); col = estilo(); r = rutas();

L = readtable(fullfile(r.dat, 'lamb_muH.csv'), 'Delimiter', ';', 'FileType', 'text');
Eexp = L.valor_meV(1); sE = L.inc_meV(1);
A = L.valor_meV(2); sA = L.inc_meV(2);
B = L.valor_meV(3); sB = L.inc_meV(3);
T = readtable(fullfile(r.dat, 'radio_proton.csv'), 'Delimiter', ';', 'FileType', 'text');
i06 = find(strcmp(T.etiqueta, 'CODATA 2006'));
r06 = T.rp_fm(i06); s06 = T.inc_fm(i06);

rp = sqrt((A - Eexp)/B);
srp = sqrt(sA^2 + sE^2 + (rp^2*sB)^2)/(2*B*rp);
rng(p.semilla);
rmc = sqrt(((A + sA*randn(p.Nmc, 1)) - (Eexp + sE*randn(p.Nmc, 1)))./(B + sB*randn(p.Nmc, 1)));

sM = escalas(c, p.asig.mL, p.asig.mN, p.asig.Z, 1);
Blo = sM.coefLO(2)*1e3;
d = corrimiento(2, 0, p.asig.Z, rp/sM.a, p.modelo, p.N, p.K);
Bnum = Blo*d.exacto/d.LO;
rpLO = sqrt((A - Eexp)/Blo);
rpNum = sqrt((A - Eexp)/Bnum);
Epred = A - B*r06^2;
sEpred = sqrt(sA^2 + (r06^2*sB)^2 + (2*B*r06*s06)^2);
nsig = (Eexp - Epred)/sqrt(sEpred^2 + sE^2);

fid = fopen(fullfile(r.res, 'extraccion_rp.txt'), 'w', 'n', 'UTF-8');
for f = [1 fid]
    fprintf(f, 'r_p (teoria completa, propagacion lineal) = %.5f +/- %.5f fm\n', rp, srp);
    fprintf(f, 'r_p (Monte Carlo, N = %d)                = %.5f +/- %.5f fm\n', p.Nmc, mean(rmc), std(rmc));
    fprintf(f, 'B teoria completa                         = %.4f meV/fm^2\n', B);
    fprintf(f, 'B orden lider (ec. 8)                     = %.4f meV/fm^2\n', Blo);
    fprintf(f, 'B numerico exacto (esfera, no relat.)     = %.4f meV/fm^2\n', Bnum);
    fprintf(f, 'r_p usando B orden lider                  = %.5f fm\n', rpLO);
    fprintf(f, 'r_p usando B numerico                     = %.5f fm\n', rpNum);
    fprintf(f, 'DeltaE_L predicho con CODATA 2006         = %.4f +/- %.4f meV\n', Epred, sEpred);
    fprintf(f, 'DeltaE_L medido                           = %.4f +/- %.4f meV\n', Eexp, sE);
    fprintf(f, 'Diferencia medido - predicho              = %.4f meV  (%.1f sigma)\n', Eexp - Epred, nsig);
end
fclose(fid);

rr = linspace(0.80, 0.90, 300);
fig = figure('Units', 'centimeters', 'Position', [2 2 24 8.5]);
subplot(1, 2, 1);
patch([r06 - s06, r06 + s06, r06 + s06, r06 - s06], [201.6 201.6 202.8 202.8], col.claro, 'EdgeColor', 'none');
hold on
h1 = plot(rr, A - B*rr.^2, '-', 'Color', col.azul);
h2 = plot(rr, A - Blo*rr.^2, '--', 'Color', col.gris);
h3 = plot(rr([1 end]), Eexp*[1 1], '-', 'Color', col.naranja, 'LineWidth', 1.2);
plot(rp, Eexp, 'o', 'Color', col.naranja, 'MarkerFaceColor', col.naranja, 'MarkerSize', 6);
plot(r06, Epred, 's', 'Color', col.gris, 'MarkerFaceColor', col.gris, 'MarkerSize', 6);
text(r06, 201.68, 'CODATA 2006', 'HorizontalAlignment', 'center', 'Color', col.gris);
set(gca, 'Layer', 'top');
xlim([0.80 0.90]); ylim([201.6 202.8]);
xlabel('r_p (fm)'); ylabel('\DeltaE_L = E(2P_{1/2}) - E(2S_{1/2})  (meV)');
legend([h1 h2 h3], {['teor' char(237) 'a completa'], ['solo t' char(233) 'rmino de ec. 8'], ['medici' char(243) 'n CREMA']}, 'Location', 'northeast');
title(['(a) Corrimiento de Lamb en ' p.asig.nombre]);
subplot(1, 2, 2);
histogram(rmc, 80, 'Normalization', 'pdf', 'FaceColor', col.azul, 'EdgeColor', 'none', 'FaceAlpha', 0.75);
hold on
xg = linspace(min(rmc), max(rmc), 300);
plot(xg, exp(-(xg - rp).^2/(2*srp^2))/(srp*sqrt(2*pi)), '-', 'Color', col.naranja);
xlim([rp - 5*srp, rp + 5*srp]);
xlabel('r_p (fm)'); ylabel('densidad de probabilidad (fm^{-1})');
title(sprintf('(b) Monte Carlo: r_p = %.5f(%.0f) fm', mean(rmc), round(1e5*std(rmc))));
legend({'10^6 muestras', ['propagaci' char(243) 'n lineal']}, 'Location', 'northeast');
guardar_figura(fig, 'fig5_extraccion_rp');
