% Dynamic analysis of a four bar linkage system, after the kinematic analysis is completed.
% Chandler Lagarde
% Collaboration with: Matthew Begnaud, Jace Delcambre, and Ronnie Kisor



%% Kinematic analysis
% Fullrun position, velocity, and acceleration analysis for any given four
% bar linkage

% This routine solves for the unknowns theta3, theta4, omega3, omega4,
% alpha2, alpha3, and alpha4

clear
clc
options = optimset('Display','off');


global a b c d t2 o2 a2 
% This allows me to set the values of these variables in functions that are
% called without having to pass them as input parameters

% This is the given value of theta2 in degrees for the graphical solution
% that provides estimates of theta3 and theta4 
t2start=0;
% Stop function to differ in where the crank will stop rotating                     
t2stop=360;
% Measured theta3 in degrees
t3est=106.6;
% Measured theta4 in degrees
t4est=131.8;

% vector of the two initial guesses of t3 and t4 in radians
tunknownguess=[t3est; t4est]*pi/180;

% Input values for global variables, these values are given
% Link2
a=2;
% Link3
b=7;
% Link4
c=9;
% Link1
d=6;
% omega2 value
o2 = 10;
% alpha2 value, doesn't need to be changed
a2 = 0;

% In case this is non-integer, the value is truncated at whole number to
% make it an integer
t2=t2start*pi/180;
n=floor(((t2stop+5)-t2start)/5);

i=1;
unknownthetas=zeros(n,2); %preallocates memory space for the matrix
Fval=zeros(n,2); 

while (1) % Will continue until the condition to break 
                  %  out of the while loop is met
[unknownthetas(i,:),Fval(i,:)]=fsolve(@thetwoeqtns,tunknownguess);
t2=t2+5*pi/180;
tunknownguess=unknownthetas(i,:);
i=i+1;
if i>n
    break
end
end

% Generates the full vector of t2 values in radians
t2 = (t2start:5:t2stop)'*pi/180;
% Generates the full vector of t2 values in degrees
t2deg = (t2start:5:t2stop)';

%Creates arrays for t3 and t4
t3 = unknownthetas(:,1);
t4 = unknownthetas(:,2);

% Creates values of t3 and t4 in degrees
t3deg = unknownthetas(:,1)*180/pi;
t4deg = unknownthetas(:,2)*180/pi;

omegas = zeros (2,n);
for i = 1:n
    A = [-b*sin(t3(i)) , c*sin(t4(i));...
        b*cos(t3(i)) , -c*cos(t4(i))];
    B = a*o2*[sin(t2(i)) , -cos(t2(i))]';
    omegas(:,i) = A\B;
end

% Creates arrays for omega3 and omega4
o3 = omegas(1,:)';
o4 = omegas(2,:)';

alphas = zeros(2,n);
for i = 1:n
    C = [-b*sin(t3(i)) , c*sin(t4(i));...
        b*cos(t3(i)) , -c*cos(t4(i))];
    D = a*a2*[sin(t2(i)) , -cos(t2(i))]' + a*o2^2*[cos(t2(i)), sin(t2(i))]' + b*(o3(i))^2*[cos(t3(i)), sin(t3(i))]' - c*(o4(i))^2*[cos(t4(i)), sin(t4(i))]';
    alphas(:,i) = C\D;
end

a3 = alphas(1,:)';
a4 = alphas(2,:)';

% % Plot
% subplot(2,2,1)
% plot(t2deg, o3) %t2 vs o3
% xlabel('theta2, deg')
% ylabel('omega3, rad/s')
% axis tight
% grid on 
% subplot (2,2,2)
% plot(t2deg, o4) %t2 vs o4
% xlabel('theta2, deg')
% ylabel('omega4, rad/s')
% axis tight
% grid on 
% subplot (2,2,3)
% plot(t2deg, a3) %t2 vs a3
% xlabel('theta2, deg')
% ylabel('alpha3, rad/s')
% axis tight
% grid on 
% subplot (2,2,4)
% plot(t2deg, a4) %t2 vs a4
% xlabel('theta2, deg')
% ylabel('alpha4, rad/s')
% axis tight
% grid on



%% Dynamic Analysis

% Create variables for the dynamic analysis
m2 = a*10^-2/0.25;
m3 = b*10^-2/0.25;
m4 = c*10^-2/0.25;
IO2 = (1/3)*m2*a^2;
IG3 = (1/12)*m3*b^2;
IO4 = (1/3)*m4*c^2;

% Acceleration of link 2 in the x and y direcitons
accelXG2 = -(a/2)*o2^2*cos(t2) - (a/2)*a2*sin(t2);
accelYG2 = -(a/2)*o2^2*sin(t2) + (a/2)*a2*cos(t2);

% Acceleration of link 3 in the x and y direcitons
accelXG3 = a*-cos(t2)*o2^2 - a*sin(t2)*a2 + (b/2)*(-cos(t3).*o3.^2 - sin(t3).*a3);
accelYG3 = a*(-sin(t2)*o2^2 + cos(t2)*a2) + (b/2)*(-sin(t3).*o3.^2 + cos(t3).*a3);

% Acceleration of link 4 in the x and y direcitons
accelXG4 = (c/2)*(-cos(t4).*o4.^2 - sin(t4).*a4);
accelYG4 = (c/2)*(-sin(t4).*o4.^2 + cos(t4).*a4);

% Dynamic Equations in matrix form
for i = 1:n
    dyMat = [1 0 0 a*sin(t2(i)) a*cos(t2(i)) 0 0 0 0 ; ...
             0 1 0 -1 0 0 0 0 0 ; ...
             0 0 1 0 1 0 0 0 0 ; ...
             0 0 0 1 0 -1 0 0 0 ; ...
             0 0 0 0 -1 0 1 0 0 ; ...
             0 0 0 (b/2)*sin(t3(i)) (b/2)*cos(t3(i)) (b/2)*sin(t3(i)) (b/2)*cos(t3(i)) 0 0 ; ...
             0 0 0 0 0 -c*cos(t4(i)) b*sin(t4(i)) 0 0 ; ...
             0 0 0 0 0 1 0 -1 0 ; ...
             0 0 0 0 0 0 -1 0 1 ];
    knownsMat = [ IO2*a2 ; m2*accelXG2(i) ; m2*accelYG2(i) ; m3*accelXG3(i) ; m3*accelYG4(i) ; IG3*a3(i) ; IO4*a4(i) ; m4*accelXG4(i) ; m4*accelYG4(i) ];
    
    dyUnknowns(:,i) = dyMat\knownsMat;
end

% Takes the dyUnknowns matrix and splits it into single variables
torque = dyUnknowns(:,1);
F1x = dyUnknowns(:,2);
F1y = dyUnknowns(:,3);
F2x = dyUnknowns(:,4);
F2y = dyUnknowns(:,5);
F3x = dyUnknowns(:,6);
F3y = dyUnknowns(:,7);
F4x = dyUnknowns(:,8);
F4y = dyUnknowns(:,9);
