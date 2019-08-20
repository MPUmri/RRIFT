function [pkParams, resid, estKepRR, stdKepRR, p] = CLRRM(Ct, Crr, t, kepRef)
  % Constrained Linear Reference Region fit
  % If kepRef > 0, then this value will be used as kepRef
  %           = 0, then mean will be used as kepRef
  %           = -1, then median will be used
  %           = -2, then interquartile mean will be used as kepRef
  
    % Use interquartile mean by default
    if nargin<4
        kepRef = -1;
    end
%%
    if kepRef <=0
        % Do first run of LRRM
        p = LRRM(Ct,Crr,t,1);
        p1 = p(:,1); % First fitted parameter
        p2 = p(:,2); % Second fitted parameter
        goodVals = (p1>0) & (p2>0) & p(:,3)>0; % Good values are when all parameters are positive
        x= p2(goodVals)./p1(goodVals); % Get kepRR, which is ratio of p2/p1
        if kepRef == 0
            % Using the mean of positive values
            estKepRR = mean(x);
            stdKepRR = std(x);
        elseif kepRef == -1
            % Using the median of positive values
            estKepRR = median(x(x>0));
            qtRange = quantile(x,[.25 .75]);
            stdKepRR = qtRange(2)-qtRange(1); % Use interquartile range as approximate for stdKepRR
        elseif kepRef == -2
            % Using the interquartile mean of the positive values
            qtRange = quantile(x,[.25 .75]);
            xMask = x>qtRange(1) & x<qtRange(2);
            estKepRR = mean(x(xMask));
            stdKepRR = std(x(xMask));
        end
    else
        estKepRR = kepRef;
        stdKepRR = nan;
        p = nan;
    end
    
    stepSize = t(2)-t(1); % Assuming constant step size throughout acquisition
    
    % Initialize matrices
    [sT, sX] = size(Ct);
    
    M1 = zeros(sT,1);
    M2 = M1;

    M1 = Crr;
    M2 = stepSize*cumtrapz(Crr);
    
    resid = zeros(sX,1);
    pkParams = zeros(sX,2);

    % Solve the linear problem
    % Disabling warnings for speed (possible non-tissue regions in image give warnings)
    warning off
    
    for i=1:sX
        curCt = squeeze(Ct(:,i));
        M3 = -stepSize*cumtrapz(curCt);
        M = [M1+estKepRR*M2, M3];
        pkParams(i,:) = mldivide(M,curCt);
        resid(i) = norm(curCt-M*pkParams(i,:)');
    end
    warning on
    
    pkParams(:,3)=pkParams(:,2);
    pkParams(:,2)=estKepRR*pkParams(:,1);
    
    % pkParams = [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kTrans_TOI/ve_TOI]
    % desire: pkParams = [kTrans_TOI/ktrans_RR, ve_TOI/ve_RR, kTrans_TOI/ve_TOI]
    pkParams(:,2) = pkParams(:,2)./pkParams(:,3);
end