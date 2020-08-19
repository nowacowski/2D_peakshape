%{
    Perform fitting of sum of multivariate normal distribution functions to
    series of 2D peaks.
    
    Provide 
    - 3D array Data: series of 2D amplitude maps of peaks
    - vectors X, Y, T for x-,y-,z-axis
    
    Use 'prepare_input.m' to make sure your data has correct format.
    Run 'Fit2D_gui.m' GUI program to perform and visualize fitting.
    
    *use of 'Fit2D.m' described better in the file itself.
%}
%% create series of 2D peaks
clear

X = linspace(14,15,50);
Y = linspace(14,15,50);

[X2, Y2] = meshgrid(X,Y);
xdata(:,:,1) = X2; 
xdata(:,:,2) = Y2;

A = -0.5;
mu_x = 14.3;
mu_y = 14.4;
sigma_x = 0.15;
sigma_y = 0.15;
rho = 0.2;

x1 = [A, mu_x, sigma_x, mu_y, sigma_y, rho];
x2 = x1 + [0.05, 0.5, 0, 0.3, 0.03, 0.1];
x = [x1, x2];

T = 1:10;
Data = zeros(length(Y), length(X), 10);
for i = 1:length(T)
    x0 = x + [0.01, 0.01, 0, 0.01, 0.01, 0, 0.01, 0.01, 0, 0.01, 0.01, 0]*i;
    dat = Gauss2D(x, xdata);
    dat_noise = awgn(dat, 20, 'measured');
    Data(:,:,i) = dat_noise;
end
%% Fit

% Prepare input for fitting
data = prepare_input(Data, X, Y, T);

% Open GUI to perform fitting
[results, errors, tw] = Fit2D_gui(data);

%% Create model from parameters obtained from fitting
model = zeros(size(Data));
for i = 1:size(Data,3)
    model(:,:,i) = Gauss2D(results, xdata);
end