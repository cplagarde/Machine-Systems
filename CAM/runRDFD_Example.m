% Run SCCA design code
% apply the scale factors h, omega, and beta
%
% This is a rise dwell fall dwell problem
% The rate of rotation is 1500 RPM
% The lift is 1" over 60 degrees followed by a dwell of 60 degrees and then
% a fall of 360 degrees
%
h=1; %units are inches
omega=2*pi*1500/60; 
betarise=pi/3;
highdwell=pi/3;
betafall=pi/3;
lowdwell=2*pi-(betarise+highdwell+betafall);
%
N=100; %this is the number of points used to map the dwells
        %This is overkill, so feel free to reduce N
        %
%The SCCA function allows you to select the type of SCCA fit and whether
%it's a rise or a fall.  The latter is actually dumb, since it's easy to
% since the fall segments can be constructed from the rise functions as I
% do in this example
%
[xrise yrise yprise ydblprise ytrplprise]=scca('modified sine','rise');
%
%You may want to open (BUT DO NOT MODIFY) the SCCA function to see how to
%specify the different functional fits, e.g. 'modified sine'
%
%SCCA generates row vectors that I paste together with the other dwell and 
% fall segments later in the code
yhighdwell=ones(1,N);
yphighdwell=zeros(1,N);
ydblphighdwell=zeros(1,N);
ytrplphighdwell=zeros(1,N);
%Since I'm using a modified trapezoid for both the rise and fall, I can
% speed up the code by not calling the scca function all over
%[xfall yfall ypfall ydblpfall ytrplpfall]=scca('modified trapezoid','fall')
xfall=xrise;
yfall=1-yrise;
ypfall=-yprise;
ydblpfall=-ydblprise;
ytrplpfall=-ytrplprise;
%
%Note that separate constructions of normalized high and low dwells, as is
%   done in this code is also dumb. This is another way in which my code is
%   not concise or efficient
ylowdwell=zeros(1,N);
yplowdwell=zeros(1,N);
ydblplowdwell=zeros(1,N)
ytrplplowdwell=zeros(1,N);
%
%The theta vectors corresponding to the rise, highdwell, fall, and low
%dwell segments are constructed below
thetarise=xrise*betarise;
thetahighdwell=betarise+[1:N]*highdwell/N;
thetafall=max(thetahighdwell)+betafall*xfall;
thetalowdwell=max(thetafall)+[1:N]*lowdwell/N;
%
% Below, the entire 360 degree set of theta, S, V, A, and J values are
% assembled into vectors
%
theta=[thetarise thetahighdwell thetafall thetalowdwell];
S=[yrise yhighdwell yfall ylowdwell]*h;
V=omega*h*[yprise/betarise yphighdwell ypfall/betafall yplowdwell];
A=omega^2*h*[ydblprise/betarise^2 ydblphighdwell ydblpfall/betafall^2 ydblplowdwell];
J=omega^3*h*[ytrplprise/betarise^3 ytrplphighdwell ytrplpfall/betafall^3 ytrplplowdwell];
%
%Figure 1 is generated by the SCCA function and shows only the normalized
%lift, velocity, and acceleration during the rise segment for the selected
%SCCA function type
%
figure(2)
subplot(2,2,1)
plot(theta*180/pi,S)
axis tight
xlabel('Theta, deg')
ylabel('Displacement, in')
grid on
subplot(2,2,2)
plot(theta*180/pi,V)
axis tight
xlabel('Theta, deg')
ylabel('Velocity, in/s')
grid on
subplot(2,2,3)
plot(theta*180/pi,A)
xlabel('Theta, deg')
ylabel('Acceleration, in/s^2')
axis tight
grid on
subplot(2,2,4)
plot(theta*180/pi,J)
axis tight
xlabel('Theta, deg')
ylabel('Jerk, in/s^3')
grid on
maxS=max(S)
maxV=max(V) 
maxA=max(A) 
maxJ=max(J)
