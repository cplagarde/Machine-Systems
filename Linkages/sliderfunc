% Script that forms a position, velocity, and acceleration analysis for an
% offset slider. 

clear
clc
global a d t2 o2
% Calls global variables that are given in the problem

a = 5;  % Length of link a in meters
d = 20;  % Length of link d in meters
t2start = 0;  % Start of where the angle of O2 pin starts
t2stop = t2start + 360;  % End of t2 rotation
t2 = (t2start:5:t2stop)*pi/180;  % Writes t2 into a radian array
o2 = 10;  % Given values of omega2 in rad/s
a2 = 0;  % Given value of alpha2 in rad/s^2
% Given variables

t4 = atan2 (d + a*sin(t2), a*cos(t2));
n = floor(((t2stop+5)-t2start)/5);
% t4 solves for the radian values of t4 and creates an array
% n creates the step size from start to stop divided by 5



% Position Analysis
c = zeros(1,n);
for i = 1:n
    c (:,i) = (d+a*sin(t2(:,i)))/(sin(t4(:,i)));
end



% Velocity Analysis
vel = zeros(2,n);
for i = 1:n
    A = [ cos(t4(i)) , -c(i)*sin(t4(i)) ; sin(t4(i)) , c(i)*cos(t4(i)) ];
    B = a*o2*[ -sin(t2(i)) ; cos(t2(i)) ];
    vel (:,i) = A\B;
end
cdot = vel (1,:);
o4 = vel (2,:);



% Acceleration Analysis
accel = zeros(2,n);
for i = 1:n
    C = [ cos(t4(i)) , -c(i)*sin(t4(i)) ; sin(t4(i)) , c(i)*cos(t4(i)) ];
    D = -a*o2^2*[cos(t2(i)) ; sin(t2(i))] + c(i)*o4(i)^2*[cos(t4(i)) ; sin(t4(i))] + 2*cdot(i)*o4(i)*[sin(t4(i)) ; -cos(t4(i))] ;
    accel (:,i) = C\D;
end
cdoubdot = accel (1,:);
a4 = accel (2,:);



% Plot Functions
t2deg = t2*(180/pi);  % Creates t2 in degrees
t4deg = t4*(180/pi);  % Creates t4 in degrees


subplot(3,2,1)
plot(t2deg, c) 
xlabel('theta2, deg')
ylabel('c, m')
axis tight
grid on 
% Plots c vs. theta2

subplot (3,2,2)
plot(t2deg, t4deg) 
xlabel('theta2, deg')
ylabel('theta4, deg')
axis tight
grid on 
% Plots theta4 vs. theta2

subplot (3,2,3)
plot(t2deg, cdot) 
xlabel('theta2, deg')
ylabel('cdot, m/s')
axis tight
grid on 
% Plots cdot vs. theta2

subplot (3,2,4)
plot(t2deg, o4)
xlabel('theta2, deg')
ylabel('omega4, rad/s')
axis tight
grid on
% Plots omega4 vs. theta2

subplot (3,2,5)
plot(t2deg, cdoubdot)
xlabel('theta2, deg')
ylabel('cdoubledot, m/s^2')
axis tight
grid on
% Plots cdoubledot vs. theta2

subplot (3,2,6)
plot(t2deg, a4)
xlabel('theta2, deg')
ylabel('alpha4, rad/s^2')
axis tight
grid on
% Plots alpha4 vs. theta2
