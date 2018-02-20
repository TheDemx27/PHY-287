clear all
clc
% Finding resonant frequency/compnent values of RLC circuit based on noise
% measurements of that circuit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define values
kb = 1.38064852*10^-23;             % Boltzmann's Constant (J/K)
tau = 80*10^-9;                     % Timestep
Fs = 1/tau;                         % Sampling Frequency
N = 10^6;                           % Number of samples
df = Fs/N;                          % Sample frequency step size
t = 0:tau:tau*(N-1);                % Time vector

%freqs = (-N/2:N/2-1)*df;           % Whole frequency vector
fpos = (0:N/2-1)*df;                % Positive frequency vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the input data
file = fopen('RLCout.dat');
x = fread(file,1000000,'float64');  % Input data stored in 'x'
fclose(file);

X = fft(x);                         % Take fft of input data        
Xpos = X(1:N/2);                    % Take only positive values into Xpos

% Take PSD, normalize then scale by 2 to accomadate for the lack of
% negative values.
PSD = 2*(real(Xpos).^2 + imag(Xpos).^2)/(length(x)*Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[val,inx] = max(PSD);               % Index of resonant frequency
f0 = inx*df;                        % Resonant frequency
w0Index = 10^5/df;                  % Index of resonant frequency

bandwidth = 11000;                  % Bandwidth estimate by graph inspection
R = mean(PSD(1:81))/(4*kb*320);     % Estimate of R
L = R/bandwidth;                    % Estimate of L
C = 1 / (L * (2*pi*fpos(8003))^2);  % Estimate of C

x0 = [R,L,C];                       % Initial guess of values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fdom = fpos(80000:120000);
PSDdom = PSD(80000:120000);
% Define our model, with parameters 1,2,3 -> R,L,C
model = @(params)(4*kb*320*params(1)*(1/sqrt(params(2)*params(3)))^4./((params(1)/params(2))^2*(2*pi*fdom).^2 + ((1/sqrt(params(2)*params(3))).^2 - (2*pi*fdom).^2).^2));

modelfull = @(params)(4*kb*320*params(1)*(1/sqrt(params(2)*params(3)))^4./((params(1)/params(2))^2*(2*pi*fpos).^2 + ((1/sqrt(params(2)*params(3))).^2 - (2*pi*fpos).^2).^2));

% Define our res function. The uncertainty is approximated by the current
% model, since the input to the circuit input is random noise and the output is
% approximated by our model, the larger the output from the model the
% greater the uncertainty in the measurement at that point.
model_res = @(params)(modelfull(params)-transpose(PSD))./modelfull(params);

% Perform least squares fit twice on the model.
[p0,~] = lsqnonlin(model_res,x0);
[p1,resnorm,~,~,~,~,jacobian] = lsqnonlin(model_res,p0);

jacobian = full(jacobian);              % Take the jacobian for uncertainty
unc = sqrt((jacobian'*jacobian)^(-1));  % Unpackage jacobian values

fprintf('\nR = %.3f +- %.3f Ohms \n', round(p1(1),3),round(unc(1,1),3))
fprintf('L = %.3f +- %.3f uH \n', round(p1(2)*10^6,3),round(unc(2,2)*10^6,3))
fprintf('C = %.3f +- %.4f pF \n', round(p1(3)*10^12,3),round(unc(3,3)*10^12,4))
fprintf('Chi^2 per degree freedom = %f\n',resnorm/(numel(x)));
% The Chi^2 per degree freedom is off due to me not being able to get
% lsqnonlin to work correctly. It isn't a perfect fit for sure, but its
% better than nothing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
figure(1)
plot(t,x)
title('Input Signal')
xlabel('Time (seconds)')
ylabel('Volts')

figure(2)
subplot(2,1,1)
plot(fpos,PSD)
title('Frequency Spectrum')
xlabel('Frequency (Hertz)')
ylabel('PSD (Volts squared per unit Frequency)')

subplot(2,1,2)
loglog(fpos,PSD,'.')
hold on;
loglog(fpos,modelfull(p1),'r')
xlabel('Frequency (Hertz)')
ylabel('PSD (Volts squared per unit Frequency)')
