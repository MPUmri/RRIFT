function ct = ToftsKety(Cp, pkParams, t)

    if isrow(Cp)
        Cp = Cp';
    end
    if isrow(t)
        t = t';
    end
    
    kTrans = pkParams(1);
    kep = pkParams(2);
    
    ct = kTrans * expConv(Cp,kep,t);
    
    if length(pkParams) == 3
    	ct = ct + pkParams(3)*Cp;
    end
end

