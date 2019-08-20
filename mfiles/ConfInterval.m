function [CI, AVG] = ConfInterval(x, intrvl)

if nargin < 2
    intrvl = [0.025 0.975];
end
%%
if length(size(x)) == 2
    [nx ny] = size(x);
    CI = zeros(nx,2);
    AVG = zeros(nx,1);

    for i=1:nx
        vals = x(i,:);
        SEM = nanstd(vals)/sqrt(length(vals));
        ts = tinv(intrvl, length(vals)-1);
        CI(i,:) = nanmean(vals) + ts*SEM;
        AVG(i) = nanmean(vals);
    end
elseif length(size(x))==3
    [nx ny, nz] = size(x);
    CI = zeros(nx,ny,2);
    AVG = zeros(nx,ny);

    for i=1:nx
        for j=1:ny
            vals = x(i,j,:);
            SEM = nanstd(vals)/sqrt(length(vals));
            ts = tinv(intrvl, length(vals)-1);
            CI(i,j,:) = nanmean(vals) + ts*SEM;
            AVG(i,j) = nanmean(vals);
        end
    end
end

