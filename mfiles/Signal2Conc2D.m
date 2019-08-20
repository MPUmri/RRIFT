function [ Ct ] = Signal2Conc2D(dataSignal, dataT1, TR, FlipAngle, r1, numBaseline)
    % dataT1 in ms
    % TR in ms
    % FlipAngle in degrees
    % r1 in 1/(mM s)
    %%

    %%  
    R10  = 1000./dataT1;
    TR = TR/1000;
    alpha = FlipAngle*pi/180;
    %%
    [nT, ~] = size(dataSignal);
    cosAlpha = cos(alpha);
    R10 = repmat(R10',[nT 1]);
    E0 = exp(-R10*TR);
    
    %%
    s = dataSignal ./ repmat(mean(dataSignal(1:numBaseline,:),1),[nT 1]); % Signal normalized by pre-contrast signal
    %%
    E = (1 - s + s .* E0 - E0 .* cosAlpha) ./ ...
        (1 - s*cosAlpha + s.*E0*cosAlpha - E0*cosAlpha);
    %%
    r1Data = real((-1/TR) * log(E));
    Ct = (r1Data - R10)/r1;
end

