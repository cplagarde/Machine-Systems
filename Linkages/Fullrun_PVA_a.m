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

t2start=0;
% This is the given value of theta2 in degrees for the graphical solution
% that provides estimates of theta3 and theta4                    
t2stop=360;
% Stop function to differ in where the crank will stop rotating  
t3est=106.6;
% Measured theta3 in degrees
t4est=131.8;
% Measured theta4 in degrees

tunknownguess=[t3est; t4est]*pi/180;
% vector of the two initial guesses of t3 and t4 in radians


a=2;
% Link2
b=7;
% Link3
c=9;
% Link4
d=6;
% Link1
o2 = 10;
% omega2 value
a2 = 0;
% alpha2 value, doesn't need to be changed

% Input values for global variables, these values are given


% The first calculation uses the specified starting value of theta2 and
% provides initial estimates of theta3 and theta4 from the graphical
% solution
%
% After that, I am going to increment theta2 by five degrees, and use the
% values of theta3 and theta4 from the previous calculation as the
% initial estimates for the new ones
% 
%  I'm going to repeat (iterate) this procedure until I've returned back to


t2=t2start*pi/180;
n=floor(((t2stop+5)-t2start)/5);
% In case this is non-integer, the value is truncated at whole number to
% make it an integer

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

t2 = (t2start:5:t2stop)'*pi/180;
% Generates the full vector of t2 values in radians
t2deg = (t2start:5:t2stop)';

% Generates the full vector of t2 values in degrees
t3 = unknownthetas(:,1);
t4 = unknownthetas(:,2);
%Creates arrays for t3 and t4

t3deg = unknownthetas(:,1)*180/pi;
t4deg = unknownthetas(:,2)*180/pi;
% Creates values of t3 and t4 in degrees





omegas = zeros (2,n);
for i = 1:n
    A = [-b*sin(t3(i)) , c*sin(t4(i));...
        b*cos(t3(i)) , -c*cos(t4(i))];
    B = a*o2*[sin(t2(i)) , -cos(t2(i))]';
    omegas(:,i) = A\B;
end

o3 = omegas(1,:)';
o4 = omegas(2,:)';
% Creates arrays for omega3 and omega4





alphas = zeros(2,n);
for i = 1:n
    C = [-b*sin(t3(i)) , c*sin(t4(i));...
        b*cos(t3(i)) , -c*cos(t4(i))];
    D = a*a2*[sin(t2(i)) , -cos(t2(i))]' + a*o2^2*[cos(t2(i)), sin(t2(i))]' + b*(o3(i))^2*[cos(t3(i)), sin(t3(i))]' - c*(o4(i))^2*[cos(t4(i)), sin(t4(i))]';
    alphas(:,i) = C\D;
end

a3 = alphas(1,:)';
a4 = alphas(2,:)';






subplot(2,2,1)
plot(t2deg, o3) %t2 vs o3
xlabel('theta2, deg')
ylabel('omega3, rad/s')
axis tight
grid on 
subplot (2,2,2)
plot(t2deg, o4) %t2 vs o4
xlabel('theta2, deg')
ylabel('omega4, rad/s')
axis tight
grid on 
subplot (2,2,3)
plot(t2deg, a3) %t2 vs a3
xlabel('theta2, deg')
ylabel('alpha3, rad/s')
axis tight
grid on 
subplot (2,2,4)
plot(t2deg, a4) %t2 vs a4
xlabel('theta2, deg')
ylabel('alpha4, rad/s')
axis tight
grid on
