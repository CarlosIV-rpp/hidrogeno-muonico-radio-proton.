function guardar_figura(fig, nombre)
r = rutas();
archivo = fullfile(r.fig, nombre);
try
    exportgraphics(fig, [archivo '.pdf'], 'ContentType', 'vector');
    exportgraphics(fig, [archivo '.png'], 'Resolution', 200);
catch
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits', 'centimeters', 'PaperSize', pos(3:4), ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 pos(3:4)]);
    print(fig, archivo, '-dpdf', '-painters');
    print(fig, archivo, '-dpng', '-r200');
end
end
