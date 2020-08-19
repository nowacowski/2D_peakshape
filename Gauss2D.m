function F = Gauss2D(x,xdata,border)

if nargin < 3
    border = ones(size(xdata(:,:,1)));
end

nx = length(x)/6;
F3 = zeros(size(xdata,1),size(xdata,2),nx);
for i = 1:nx
    mux = x(6*i-4); muy = x(6*i-2);
    sigx = x(6*i-3); sigy = x(6*i-1);
    rho = x(6*i);
    amp = x(6*i-5);
    factor = -1/(2*(1-rho^2));
    term1 = ((xdata(:,:,1)-mux)/sigx).^2;
    term2 = ((xdata(:,:,2)-muy)/sigy).^2;
    term3 = (2*rho*(xdata(:,:,1)-mux).*(xdata(:,:,2)-muy))/(sigx*sigy);
    F3(:,:,i) = amp*exp(factor*(term1+term2-term3));
end
F2 = sum(F3,3);

F = F2.*border;