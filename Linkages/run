% Finds values close to theta3 and theta4

clear
clc
% Example 4-1 pg 190 of Norton, 5th ed
t3est=90*pi/180;
t4est=118*pi/180;
theta3_4_guess=[t3est; t4est];
global a b c d t2
t2=-70*pi/180;
a=10;
b=10;
c=10;
d=20;

[theta3_4,Fval]=fsolve('thetwoeqtns',theta3_4_guess)
% Print calculated values to command screen in degrees
%  represent the answer to two decimal places
round(100*theta3_4'*180/pi)/100
% Note that I've taken the transpose of the column vector
%  so that the two values are appear side by side in a row 
