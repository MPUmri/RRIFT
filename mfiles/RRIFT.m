function [estKtRR, num, denum] = RRIFT(CpTail, CrrTail, tTail, kepRR)
%%
    num = CrrTail - CrrTail(1) + kepRR * cumtrapz(tTail,CrrTail);
    denum = cumtrapz(tTail,CpTail);
    if isrow(num)
       num=num'; 
    end
    if isrow(denum)
        denum=denum';
    end
    estKtRR = denum\num;
end

