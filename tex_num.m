function s = tex_num(x, cifras)
if x == 0
    s = '$0$';
    return
end
e = floor(log10(abs(x)));
if e >= -2 && e <= 3
    s = sprintf(['%.' num2str(max(cifras - 1 - e, 0)) 'f'], x);
else
    m = round(x/10^e, cifras - 1);
    if abs(m) >= 10
        m = m/10;
        e = e + 1;
    end
    s = sprintf(['%.' num2str(cifras - 1) 'f\\times10^{%d}'], m, e);
end
s = ['$' s '$'];
end
