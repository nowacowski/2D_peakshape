function [cov_err] = cov_error(resnom,xData,x,J)
ssq = resnom/(length(xData)-length(x));     %"mean residual variance"
cov2 = inv(J'*J)*ssq;                         %"variance-covariance matrix"
cov_err = full(cov2);
end