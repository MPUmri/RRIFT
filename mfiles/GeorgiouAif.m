function [C_p, C_b] = GeorgiouAif(t, t0)
% Quick script for generating a population-based AIF...

    if nargin < 1
        t=0:499;
        t=t'/60;
    end
    if nargin < 2
        t0 = 0;
    end
    if isrow(t)
        t = t';
    end
%%
    t = t-t0;
    t(t<0) = 0;

    % %% Population-averaged parameters from Parker
    A = [0.37, 0.33, 10.06];
    m = [0.11, 1.17, 16.02];
    alpha = 5.26;
    beta = 0.032;
    tau = 0.129;

    nList = [];
    C_b1 = zeros(length(t),1); % Exponential part
    C_b2 = C_b1;
    
    for z=1:length(t)
        C_b1(z) = sum(A.*exp(-m*t(z)));
        % Gamma variate part can't be computed for high N (when t>7.5 min)
        % so fix it to a single value
        % The value was determined by plotting C_b2 vs t
        if t(z)>7.5
            C_b2(z)=3.035;
            continue;
        end
        N = floor(t(z)/tau);
        for j=0:N
            curGamma = GammaVariate((j+1)*alpha+j, beta, t(z)-j*tau);
            if ~isnan(curGamma)
                % NaNs are returned for large, so assume they're zero
                % (They are caused by a 0*Inf term)
                C_b2(z) = C_b2(z) + curGamma;
            end
        end
    end
    
    C_b = C_b1 .* C_b2;
    C_b(t<=0)=0;    

%     C_b is the concentration in the whole blood
%     Most models use the concentration in the blood plasma C_p
%     C_b can be converted to C_p by using hematocrit.
     Hct = 0.35; % The chosen patient had Hct of 0.35
     C_p = C_b / (1-Hct);
end
