% check and prepare input data

function data = prepare_input(Data, X, Y, T)

% check if provided input arguments are correct
if nargin < 3
    error('You need to provide at least Data, X and Y arrays');
end
if length(size(Data)) < 2
    error('Data needs to be 2- or 3-dimensional array')
end    
if size(Data,1) ~= length(Y)
    error('Size of Y needs to match 1st dimention of Data')
end
if size(Data,2) ~= length(X)
    error('Size of X needs to match 2nd dimention of Data')
end
if nargin < 4
    T = 0;
end
if size(Data,3) ~= length(T)
    error('Size of T needs to match 3rd dimention of Data')
end

% sort axis in ascendent order
if size(X,1) == 1
    [X, iX] = sortrows(X');
    X = X';
elseif size(X,2) == 1
    [X, iX] = sortrows(X);
    X = X';
end
if size(Y,1) == 1
    [Y, iY] = sortrows(Y');
elseif size(Y,2) == 1
    [Y, iY] = sortrows(Y);
end
Data = double(real(Data(iY,iX,:)));

% contruct structure array with prepared input arguments
data.Data = Data;
data.X = X;
data.Y = Y;
data.T = T;