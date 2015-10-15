%This routine solves for the unknown thetas, theta3 and theta4
% for a GCRR for a full range rotation of the crank
clear
clc
global a b c d t2 %This allows me to set the values of these
                                % variables in functions that I call
                                % without having to pass them as input
                                % parameters
t2start=-70; %This is the value of theta2 in degrees for which I have 
                      % a graphical solution that provides estimates of the
                      % corresponding theta3 and theta4

t2stop=70;  %Stop function to differ in where the crank will stop rotating                      
t3est=47;  %Measured theta3 in degrees
t4est=192; %Measured theta4 in degrees
tunknownguess=[t3est; t4est]*pi/180; % vector of the two initial guesses of 
                    % theta3 and theta4,in radians
a=10;
b=10;
c=10;
d=20;

% The first calculation uses the specified starting value of theta2 and
% provides initial estimates of theta3 and theta4 from the graphical
% solution
%
% After that, I am going to increment theta2 by five degrees, and use the
% values of theta3 and theta4 from the previous calculation as the
% initial estimates for the new ones
% 
%  I'm going to repeat (iterate) this procedure until I've returned back to
t2=t2start*pi/180
n=floor(((t2stop+5)-t2start)/5); %in case this is non-integer, the value is truncated at 
                               % whole number to make it an integer
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
t2=[t2start:5:t2stop]'*pi/180; %Generates the full vector of 
                        %  theta2 values.  This could've been done along
                        %  the way, but it would've required that the
                        %  current value of t2 be passed as an additional
                        %  parameter, instead of relying on making it a
                        %  global variable
%Plot theta3 versus theta2 and theta4 versus theta2 in separate plots
subplot(2,1,1)
plot(t2*180/pi,unknownthetas(:,1)*180/pi)
axis tight
grid on
xlabel('Theta2, deg')
ylabel('Theta3, deg')
subplot(2,1,2)
plot(t2*180/pi,unknownthetas(:,2)*180/pi)
axis tight
grid on
xlabel('Theta2, deg')
ylabel('Theta4,deg')
