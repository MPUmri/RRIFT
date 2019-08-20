function [ outData ] = unravel( inData )
% Meant to convert 4D DCE data to 2D matrix
% I got tired of writing the same command over and over, so made it a function

    [sX sY sZ sT] = size(inData);
    outData = reshape(inData,[sX*sY*sZ sT]);
    
end

