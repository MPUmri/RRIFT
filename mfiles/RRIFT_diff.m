function [estKtRR, num, denum] = RRIFT_diff(CpTail, CrrTail, tTail, kepRR)
%% Differential form of RRIFT, i.e. avoids extra integrals
    num = gradient(CrrTail,tTail(2)-tTail(1)) + kepRR * CrrTail;
    denum = CpTail;
    if isrow(num)
       num=num'; 
    end
    if isrow(denum)
        denum=denum';
    end
    estKtRR = denum\num;
end

