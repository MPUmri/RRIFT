function [fittedCRR, exitflag] = DoEMM(cRR, t)
    % cRR = measured concentration curve of reference tissue
    % t = frame times

    stepSize = t(2)-t(1);

    aifGuess = [0.01,0.01,0.01,0.01,0.01,0.01];
    
    options=optimset('Algorithm','levenberg-marquardt','display','off');
    [nParams, res, resnorm, exitflag]=lsqnonlin(@(x) cRR - real(EMM(x, t)),aifGuess,[],[],options);

    fittedCRR = real(EMM(nParams,t));
end
