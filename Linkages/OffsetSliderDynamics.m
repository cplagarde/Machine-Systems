%% Chandler Lagarde
% MCHE 363: Kinematics and Dynamics of Machine Systems
% Solves a kinematic and dynamic analysis for an offset slider with a
% horizontal block at the end. 


%% Kinematic Analysis

clear
clc

a = 0.03; %% meters
b = 0.13; %% meters
o2 = 2000*2*pi/60; %% rads/s
m2 = 2.5; %% kg
m3 = 1; %% kg
mb = 0.5; %% kg
r2 = 0.01; %% meters
r3 = b/2;
mu = 0.1; 
c = 0; %% Vertical offset
a2 = 0; %% alpha 2

t2Start = 0;
t2Stop = t2Start + 360;
n = floor(((t2Stop+5)-t2Start)/5);
t2 = (t2Start:5:t2Stop).*pi/180;

% Solves for theta3 and d
t3 = zeros(1,n);
d = zeros(1,n);
for i = 1:n
    t3(1,i) = asin((a*sin(t2(i))-c)/b);
    d(1,i) = a*cos(t2(i))-b*cos(t3(i));
end

% Solves for d-dot and omega3
C = zeros(2,n);
for i = 1:n
    A = [1 -b*sin(t3(i)) ; 0 b*cos(t3(i))];
    B = [-a*o2*sin(t2(i)) ; a*o2*cos(t2(i))];
    C(:,i) = A\B;
end
ddot = C(1,:);
o3 = C(2,:);

% Solves for d-double-dot and alpha3
D = zeros(2,n);
for i = 1:n
    A = [1 -b*sin(t3(i)) ; 0 b*cos(t3(i))];
    B = a*a2*[-sin(t2(i)) ; cos(t2(i))] - a*o2^2*[cos(t2(i)) ; sin(t2(i))] + b*o3(i).^2*[cos(t3(i)) ; sin(t3(i))];
    D(:,i) = A\B;
end
ddoubdot = D(1,:);
a3 = D(2,:);



%% Dynamic Analysis

% Inertia calculated with radius of gyration
IG3 = (m3*r3^2);

XG2 = zeros(1,n);
XG3 = zeros(1,n);
YG2 = zeros(1,n);
YG3 = zeros(1,n);
unknowns = zeros(8,n);

% For loop solves for x-double-dot-G2, y-double-dot-G2,...
% x-double-dot-G3, and y-double-dot-G3. 
% Also solves dynamic matrix equations for all 8 unknowns.
for i = 1:n
    XG2(1,i)=r2*((-sin(t2(i))*a2)-(o2^2*cos(t2(i))));
    YG2(1,i)=r2*((a2*cos(t2(i)))-(o2^2*sin(t2(i))));
    XG3(1,i)=ddoubdot(i)-((b-r3)*a3(i)*sin(t3(i)))-(b-r3)*(o3(i))^2*cos(t3(i));
    YG3(1,i)=(b-r3)*(a3(i)*cos(t3(i))-(o3(i))^2*sin(t3(i)));
    
    % if loop to generate the signddot with the first value replaced with
    % zero, since the first value will be null otherwise. 
    if i == 1
    signddot(i) = 0;
    else 
    signddot(:,i) = ddot(i)/abs(ddot(i));
    end
    
    LAns = [1 0 0 a*sin(t2(i)) -a*cos(t2(i)) 0 0 0 ; ...
            0 1 0 -1 0 0 0 0 ; ...
            0 0 1 0 -1 0 0 0 ; ...
            0 0 0 1 0 -1 0 0 ; ...
            0 0 0 0 1 0 1 0 ; ...
            0 0 0 -r3*sin(180-t3(i)) -r3*cos(180-t3(i)) -(b-r3)*sin(180-t3(i)) (b-r3)*cos(180-t3(i)) 0 ; ...
            0 0 0 0 0 1 0 -mu*signddot(i) ; ...
            0 0 0 0 0 0 -1 1 ];
    RAns = [0; m2*XG2(i); m2*YG2(i); m3*XG3(i); m3*YG3(i); IG3*a3(i); mb*ddoubdot(i); mb*0];
    unknowns(:,i) = LAns\RAns;
end

% Assigns the unknowns to their respective unknown values.
Tau = unknowns(1,:);
F1X = unknowns(2,:);
F1Y = unknowns(3,:);
F2X = unknowns(4,:);
F2Y = unknowns(5,:);
F3X = unknowns(6,:);
F3Y = unknowns(7,:);
N = unknowns(8,:);

% Magnitude of the forces.
F1 = sqrt(F1X.^2+F1Y.^2);
F2 = sqrt(F2X.^2+F3Y.^2);
F3 = sqrt(F3X.^2+F3Y.^2);



%% Plotting

t2deg = t2*(180/pi);  % Creates t2 in degrees

subplot (2,2,1)
plot(t2deg, F1)
xlabel('theta2, deg')
ylabel('Force 1, N')
axis tight
grid on

subplot (2,2,2)
plot(t2deg, F2)
xlabel('theta2, deg')
ylabel('Force 2, N')
axis tight
grid on

subplot (2,2,3)
plot(t2deg, F3)
xlabel('theta2, deg')
ylabel('Force 3, N')
axis tight
grid on

subplot (2,2,4)
plot(t2deg, Tau)
xlabel('theta2, deg')
ylabel('Tau, N*m')
axis tight
grid on
