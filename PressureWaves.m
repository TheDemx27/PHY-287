% Pressure waves in a 1-d "tube".
clear all
clc

Length = 1;            % length of rod in meters

c = 345;  %speed of waves in atmosphere
lambda = Length;   %set wavelength
f = c/lambda;        %set frequency

%Spacial step size of 0.01 meters
dx = 0.01;

% Question 2

% A: The smaller the alpha the smaller the step size
% in time. alpha = 0.25 makes for a smaller step size than alpha = 1. alpha
% = 1.0025 makes almost all the values be 0.

% B: It appears that any value in the range (0,1] will give reliable
% numerical solutions.

% C: It becomes an issue when alpha >= 1 because by rearranging the
% relation, we have that deltaT = deltaX*sqrt(alpha)/c. With alpha = 1, deltaT
% will be the time it takes for the wave to go distance deltaX. If the wave
% travels a distance greater than deltaX in time deltaT, then the values
% needed to calculate the next point get shifted so that u(j,n) =
% u(j+1,n+1), which when the initial values are set to zero can render the
% whole solution to just be zero.

alpha = input('Input value for alpha: ');

dt = sqrt(alpha*dx^2/c^2);
Total_time = 8/f;       % total time for simulation (in sec)

X = 0:dx:Length;       % space grid
T = 0:dt:Total_time;   % time grid

M = length(X);         % number of grid points in length
N = length(T);         % number of steps in time

u=zeros(M,N);

u(:,1)=0;    % set initial pressure (uniform) at t=0
u(M,:)=0;     % set left boundary at atmospheric pressure in Pa

for n=2:N       %time loop
    u(1,n) = 0.001*sin(2*pi*f*n*dt); %Forcing function
    for j=2:M-1     %space loop
        if n == 2
            u(j,n) = (1-alpha)*u(j,1)+alpha/2*(u(j+1,1)+u(j-1,1));
        else
            u(j,n) = 2*(1-alpha)*u(j,n-1)-u(j,n-2)+alpha*(u(j+1,n-1)+u(j-1,n-1));
        end
    end
    %u(M,n) = u(M-1,n);  % implements Neumann boundary condition on right side.
end

figure(1)
mesh(T,X,u)
xlabel('Time in seconds');
ylabel('Position in Meters');
zlabel('Pressure Difference in Pascals');
