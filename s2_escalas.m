clear; close all; clc;
c = constantes(); p = parametros(); r = rutas();

atomos = {p.ref, p.asig};
V = zeros(14, 2);
for k = 1:2
    A = atomos{k};
    s = escalas(c, A.mL, A.mN, A.Z, p.rp);
    psiNum = zeros(2, 1);
    Pin = zeros(2, 1);
    for n = 1:2
        g = malla(1e-4*s.rhoP, 20*n/A.Z, p.N);
        v0 = potencial(g.rho, A.Z, 0, 'puntual');
        [~, u] = autovalor(n, 0, A.Z, g, v0);
        psiNum(n) = (u(1)/g.rho(1))^2/(4*pi)/s.a^3;
        dentro = g.rho <= s.rhoP;
        Pin(n) = trapz(g.x(dentro), u(dentro).^2.*g.rho(dentro));
    end
    dS = corrimiento(2, 0, A.Z, s.rhoP, p.modelo, p.N, p.K);
    dP = corrimiento(2, 1, A.Z, s.rhoP, p.modelo, p.N, p.K);
    V(:, k) = [s.mr_me; s.a; s.En(1); s.En(2); s.psi0(1); s.psi0(2); psiNum(2); ...
        s.rhoP; Pin(2); s.coefLO(2); dS.LO*s.Eh; dS.exacto*s.Eh; ...
        abs(dS.exacto*s.Eh/s.En(2)); dP.exacto*s.Eh];
end
cociente = V(:, 2)./V(:, 1);

nombres = {'m_r/m_e'; 'a (fm)'; 'E_1 (eV)'; 'E_2 (eV)'; '|psi_100(0)|^2 (fm^-3)'; ...
    '|psi_200(0)|^2 analitico (fm^-3)'; '|psi_200(0)|^2 numerico (fm^-3)'; 'r_p/a'; ...
    'P(r<r_p) en 2S'; 'coef. LO de r_p^2 en 2S (eV/fm^2)'; 'DeltaE_2S LO (eV)'; ...
    'DeltaE_2S exacto (eV)'; 'DeltaE_2S/|E_2|'; 'DeltaE_2P exacto (eV)'};
tabla = table(nombres, V(:, 1), V(:, 2), cociente, ...
    'VariableNames', {'cantidad', p.ref.nombre, 'asignado', 'cociente'});
writetable(tabla, fullfile(r.res, 'tabla_escalas.csv'), 'Delimiter', ';');
disp(tabla)

etiquetas = {'$m_r/m_e$'; '$a$ (fm)'; '$E_1$ (eV)'; '$E_2$ (eV)'; ...
    '$|\psi_{100}(0)|^2$ (fm$^{-3}$)'; '$|\psi_{200}(0)|^2$ anal\''itico (fm$^{-3}$)'; ...
    '$|\psi_{200}(0)|^2$ num\''erico (fm$^{-3}$)'; '$r_p/a$'; '$P(r<r_p)$ en 2S'; ...
    'Coeficiente de $r_p^2$ en 2S (eV/fm$^2$)'; '$\Delta E_{2S}$ orden l\''ider (eV)'; ...
    '$\Delta E_{2S}$ num\''erico exacto (eV)'; '$\Delta E_{2S}/|E_2|$'; '$\Delta E_{2P}$ num\''erico (eV)'};
fid = fopen(fullfile(r.res, 'tabla_escalas.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', '\begin{tabular}{lccc}', '\toprule', ...
    ['Cantidad & ' p.ref.tex ' & ' p.asig.tex ' & ' p.asig.tex '/' p.ref.tex ' \\'], '\midrule');
for k = 1:numel(etiquetas)
    fprintf(fid, '%s & %s & %s & %s \\\\\n', etiquetas{k}, tex_num(V(k, 1), 4), ...
        tex_num(V(k, 2), 4), tex_num(cociente(k), 4));
end
fprintf(fid, '%s\n', '\bottomrule', '\end{tabular}');
fclose(fid);
