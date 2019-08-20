function ct = EMM(aifParams, t)
    A = aifParams(1);
    q = aifParams(2);
    alpha = aifParams(3);
    beta = aifParams(4);
    gamma = aifParams(5);
    t0 = aifParams(6);
    t=t-t0;
    t(t<0)=0;

    ct = A * (1-exp(-alpha*(t))).^q ...
            .* exp(-beta*(t)) ...
            .* (1+exp(-gamma*(t)))/2;

end
