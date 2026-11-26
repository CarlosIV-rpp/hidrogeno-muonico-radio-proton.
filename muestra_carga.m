function X = muestra_carga(modelo, N)
switch lower(modelo)
    case 'esfera'
        r = sqrt(5/3)*rand(N, 1).^(1/3);
    case 'exponencial'
        r = -sum(log(rand(N, 3)), 2)/sqrt(12);
    case 'gauss'
        X = randn(N, 3)/sqrt(3);
        return
    otherwise
        error('Modelo de carga desconocido: %s', modelo);
end
u = randn(N, 3);
u = bsxfun(@rdivide, u, sqrt(sum(u.^2, 2)));
X = bsxfun(@times, r, u);
end
