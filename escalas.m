function s = escalas(c, mL, mN, Z, rp)
s.mr = mL*mN/(mL + mN);
s.mr_me = s.mr/c.me;
s.a = c.hbarc/(s.mr*c.alfa);
s.Eh = s.mr*1e6*c.alfa^2;
n = (1:2)';
s.En = -Z^2*s.Eh./(2*n.^2);
s.psi0 = Z^3./(pi*n.^3*s.a^3);
s.coefLO = (2/3)*Z^4*s.Eh./(n.^3*s.a^2);
s.rhoP = rp/s.a;
s.dE2S_LO = s.coefLO(2)*rp^2;
end
