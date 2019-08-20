function [x] = QuantileFilter(x, qtRange)

qt = quantile(x(:),qtRange);

x(x>qt(2)) = NaN;
x(x<qt(1)) = NaN;

end

