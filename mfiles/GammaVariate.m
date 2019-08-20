function [out] = GammaVariate(alpha, beta, t)

if t<0
    out = 0;
    return;
end

num = t.^alpha .* exp(-t./beta);
denum = beta^(alpha+1) * gamma(alpha+1);
out = num/denum;

end

