function [A,mu,sig] = caruana(x,y)
%Caruana's algorithm: numerical fitting of gaussian function. it linearizes 
%y = Aexp(-(x-mu)^2/2sig^2) to ln(y) = a + bx + cx^2 which is a parabola

% inputs need to be 1xn arrays
if size(x,1) > size(x,2)
    x = x';
end
if size(y,1) > size(y,2)
    y = y';
end

z = real(log(y)/[x.^2;x;ones(size(y))]);
mu = -z(2)/(2*z(1));
sig = sqrt(-1/(2*z(1)));
A = exp(z(3)-((z(2)^2)/(4*z(1))));