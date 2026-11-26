function v = potencial(rho, Z, rhoP, modelo)
if rhoP <= 0
    modelo = 'puntual';
end
switch lower(modelo)
    case 'puntual'
        v = -Z./rho;
    case 'esfera'
        R = sqrt(5/3)*rhoP;
        v = -Z./rho;
        m = rho < R;
        v(m) = -Z*(3 - (rho(m)/R).^2)/(2*R);
    case 'exponencial'
        lam = sqrt(12)/rhoP;
        y = lam*rho;
        v = -(Z./rho).*(-expm1(-y) - 0.5*y.*exp(-y));
    case 'gauss'
        sig = rhoP/sqrt(3);
        v = -Z*erf(rho/(sqrt(2)*sig))./rho;
    otherwise
        error('Modelo de carga desconocido: %s', modelo);
end
end
