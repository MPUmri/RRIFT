function [a, b] = mylls(x,y)
% function [a, b] = mylls(x,y)
%
%   a function for fast linear least squares regression of y = ones(size(x))*a + x*b
%
%   inputs:
%   x = x data arranged by column
%   y = y data arranged by column
%
%   outputs:
%   a = vector of intercepts
%   b = vector of slopes
%   

%
%   Ives Levesque, August 2012
%

n = size(x,1);

if size(y) ~= size(x)
    error('Size of x and y must match.')
end

% estimate slope
num = sum(x.*y) - sum(x).*sum(y)/n;
denom = sum(x.^2) - (sum(x).^2)/n;
b = num./denom;

% estimate intercept term
a = mean(y) - b.*mean(x);