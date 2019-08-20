function [T1, M0, E1] = doDESPOT1(flipData, alpha, TR, mask)
    % Simple DESPOT1 function - doesn't use structs
    % Input: flipData - nFlip-by-nVox data
    %        alpha - flips angle in radians
    %        TR - repetition time; also determines units of T1
    %%
    [sX sY sZ numFlips] = size(flipData);
    flipData = reshape(flipData,[sX*sY*sZ numFlips])';
    flipData = flipData(:,mask(:));
    %%
    [numFlips,nVox] = size(flipData);
    if isrow(alpha)
        alpha = alpha';
    end
    alpha = repmat(alpha,[1 nVox]);
    %%
    y = flipData ./ sin(alpha);
    x = flipData ./ tan(alpha);

    %%
    [a,b] = mylls(x,y);
    %%
    [E1, M0, T1] = deal(zeros(sX,sY,sZ));
    %%
    E1(mask(:)) = b;
    M0(mask(:)) = a ./ (1-b);
    T1(mask(:)) = real(-TR./log(E1(mask(:))));

end

