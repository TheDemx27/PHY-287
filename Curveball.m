clear all, clc

% DEFINITIONS
dt = 0.01;
tfinal = 10;
N = tfinal/dt;
t = 0:dt:tfinal;

p = zeros(3,N);
v = zeros(3,N);
w = [0,0,10*2*pi];

% USER INPUT
% theta = input('Theta = ');
% phi = input('Phi = ');

theta = 81;
phi = 16;

fprintf('theta = %d degrees \nphi   = %d degrees', theta, phi)

% Question 2: Through trial and error, 4 tuples of angles that put the ball in the goal
% are (81,16),(81,17),(81,18),(80,19)

% Question 3a: The path of the ball would be a parabola of the form
% z(t) = -1/2*9.8*t^2 + 30*sin(phi*pi/180)*t, 
% x(t) = 30*cos(theta*pi/180)*cos(phi*pi/180)*t
% y(t) = 30*sin(theta*pi/180)*cos(phi*pi/180)*t

% Quesetion 3c: The path of the ball is affected by lift, which depends of
% the relative directions and magnitudes of the velcocity and spin of the
% ball. The side of the ball that moves with the oncoming wind has a lower
% pressure on that side, the side of the ball that moves against the
% oncoming wind has a higher pressure. If we let v represent the velocity
% of the ball and w the spin vector, then v x w always points in the
% direction of the lift generated by the spin, and this always coincides
% with our intuition.

% INITIAL CONDITIONS
v(1,1) = 30*cos(theta*pi/180)*cos(phi*pi/180);
v(2,1) = 30*sin(theta*pi/180)*cos(phi*pi/180);
v(3,1) = 30*sin(phi*pi/180);

% MAIN LOOP
for ii = 1:N
    A = CalcForce(w,v(:,ii));
    vmid = v(:,ii) + A*dt/2;
    Amid = CalcForce(w,vmid);
    
    v(:,ii+1) = v(:,ii) + Amid*dt;
    p(:,ii+1) = p(:,ii) + v(:,ii)*dt;
    
    if p(1,ii+1) <= 0 || p(2,ii+1) <= 0 || p(3,ii+1) <= 0
        outcome = CheckForOutcomes(p(:,ii+1));
        break
    end
end

outcomes = ["entered the goal without touching the ground." "hit the ground without entering the goal." "missed the goal."];

% Question 3b
stem3(p(1,:),p(2,:),p(3,:))
fprintf('\nThe ball %s\n',outcomes(outcome))
axis([0 3 0 50 0 7])

function result = CheckForOutcomes(p)
    ballr = 0.10981691073;
    result = 3; %Missed the goal
    if(p(3) <= 0)
        result = 2; %Hit ground without getting into goal
    end
    if(p(1) <= 0) || p(2) <= 0
        if p(2) > 30.44 + ballr && p(2) < 37.66 - ballr
            if p(3) > 0 + ballr && p(3) < 2.66 - ballr
                result = 1; %Made the goal
            end
        end
    end
end

function A = CalcForce(w,v)
    rho = 1.22;
    Cpd = 0.3;
    g = 9.8;
    mball = 0.43;
    rball = 0.1098;
    mda = rho*4*pi/3*rball^3;
    Fl = 3/(4*pi)*mda*cross(w,v);
    Fd = -1/2*rho*Cpd*(pi*rball^2)*norm(v)*v';
    Fg = [0 0 -mball*g];
    
    A = ((Fg + Fd + Fl)/mball)';
end