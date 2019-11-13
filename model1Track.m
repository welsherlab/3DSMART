%% Simulate Tracking Performance
%This code simulates 2D real-time tracking using experimental piezo response
%% Set simulation parameters
tau = 20e-6;            %Time step [sec]
N = 1e4;                %Number of trajecotry steps
time = tau*(1:N);       %create a time vector for plotting

%% Laser scan parameters
w = .2e-6;              %Beam width
bigW = [w^2 0;0 w^2];   %2D width of the PSF
beamAmp = 0.5e-6;       %Laser Scan Size

%% Feedback Parameters
K_P = 0;                %Proportional Feedback Gain
K_I = 0.012;            %Integral Feedback Gain

%% Set diffusion  and particle parameters
D = 4e-12;              %Diffusion coefficient [m^2/sec]
k = sqrt(2*D*tau);

brightness = 6e4;       %Particle brightnesss [Hz]
bg = 700;               %Background [Hz]
%% Generate Trajectory
%Generate steps in x (dx) and y (dy)
dx = k * randn(N,1);
dy = k * randn(N,1);

%Generate XY trajectory
x = cumsum(dx);
y = cumsum(dy);
%% Initialize Position Estimation and Stage Variables

%Initial calculated position
x_k_1 = [0;0];

%Initial uncertainty
sigma_k_1 = (0.2e-6)^2;      

%Initial stage positions
prevStagePosX = 0;
prevStagePosY = 0;

%Initialize variables
err = [];
stageCommands = [];
x_k = [0 0];
out = [];
%% Load impulse response data
measuredImpulseResponse = dlmread('Impulse_Response.txt');
%% Track
% Loop through trajectory steps
for i=1:length(x)
    
    %Retrieve current beam positions
    [beamX(i),beamY(i)] = kt(i,beamAmp);
    
    %Correct current laser position with stage
    laserPositionX = beamX(i) + prevStagePosX;
    laserPositionY = beamY(i) + prevStagePosY;
    laserPosition = [laserPositionX laserPositionY]';
    
    %Get current particle positions
    particlePositionX = x(i);
    particlePositionY = y(i);
    particlePosition = [particlePositionX particlePositionY]';
    
    %Get number of photons emitted at current step
    emission(i) = poissrnd(tau*emRate(particlePosition-[prevStagePosX prevStagePosY]',laserPosition-[prevStagePosX prevStagePosY]',bg,brightness,bigW));
    
    %Update step
    x_k = (w^2*x_k_1+emission(i)*sigma_k_1*(laserPosition-[prevStagePosX prevStagePosY]'))/(w^2+emission(i)*sigma_k_1);
    sigma_k = w^2*sigma_k_1/(w^2+emission(i)*sigma_k_1);
    
    %Error calculation.
    err(i,:) = x_k;
    pkOut(i) = sigma_k;
    
    %Prediction step
    x_k_1 = x_k;
    sigma_k_1 = sigma_k + 2*D*tau;
    
    %Store position estimates
    out(i,:) = x_k;
    
    %Get current stage commands and positions
    stageCom = K_P * err(i,:) + K_I * sum(err,1);
    stageCommands(i,:) = stageCom;
    stagePosX = impulseResponse(stageCommands(:,1),measuredImpulseResponse);
    stagePosY = impulseResponse(stageCommands(:,2),measuredImpulseResponse);
    
    prevStagePosX = stagePosX(end);
    prevStagePosY = stagePosY(end);
end

%% Trajectory Statistics
% Calculate the total tracking error
%%
trackingError = (sum((x-(stagePosX+out(:,1))).^2)+sum((y-(stagePosY+out(:,2))).^2));

%Calculate emission rate
rate = (sum(emission)/max(time));

%This is the total tracking error per Hz
normalizedTrackingError = (sum((x-stagePosX).^2)+sum((y-stagePosY).^2))/rate;

%% Plot Results
figure(1)
subplot(3,1,1)
plot(time,x,time,stagePosX+out(:,1),time,stagePosX)
legend('Particle Position','Calculated Position','Stage Position')
xlabel('Time (s)')
ylabel('X Position (m)')

subplot(3,1,2)
plot(time,y,time,stagePosY+out(:,2),time,stagePosY)
legend('Particle Position','Calculated Position','Stage Position')
xlabel('Time (s)')
ylabel('Y Position (m)')
subplot(3,1,3)
plot(time,emission)
xlabel('Time (s)')
ylabel('Intensity (cts/bin)')

%% Impulse Response Function
function out = impulseResponse(x,measuredImpulseResponse)
%Convolve with stage command data
temp = conv(x,measuredImpulseResponse);
%Crop data
out = temp(1:length(x));
end

%% Knight's Tour Coordinates
function [xOut,yOut]=kt(index,d)
x = d*[-1;-0.5;-1;0;1;0.5;1;0;0.5;1;0;-1;-0.5;0.5;1;0.5;-0.5;-1;0;1;0.5;-0.5;0;-0.5;-1];
y = d*[1;0;-1;-0.5;-1;0;1;0.5;-0.5;0.5;1;0.5;-0.5;-1;0;1;0.5;-0.5;-1;-0.5;0.5;1;0;-1;0];
xOut = x(mod(index-1,length(x))+1);
yOut = y(mod(index-1,length(x))+1);
end

%% Emission
function out=emRate(x,c,b,s,w)
%x - position of particle
%c - position of beam
%b - background
%w - spatial covariance
out=s*exp(-0.5*(x-c)'*inv(w)*(x-c))+b;
end