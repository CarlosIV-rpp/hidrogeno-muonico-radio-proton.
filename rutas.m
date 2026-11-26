function r = rutas()
r.base = fileparts(fileparts(mfilename('fullpath')));
r.fig = fullfile(r.base, 'figuras');
r.res = fullfile(r.base, 'resultados');
r.dat = fullfile(r.base, 'datos');
if ~exist(r.fig, 'dir'), mkdir(r.fig); end
if ~exist(r.res, 'dir'), mkdir(r.res); end
end
