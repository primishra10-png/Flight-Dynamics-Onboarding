%% 3-DOF Rocket Simulation 
clear; clc; close all;

%% Given Parameters
thrust = 2000; % N
mass = 40; % kg
burnTime = 5; % sec
diameter = 0.1; % m
dragCoeff = 0.75; % dimensionless
g = 9.81; % m/s^2
rho = 1.225; % kg/m^3 (air density)

%% Simulation Setup
dt = 0.01; % time step (seconds)
maxTime = 60; % max simulation time
t = 0:dt:maxTime; % time array
n = length(t); % number of time steps

%% Initial Conditions
launchAngle = 85; % degrees (this is an assumption because no angle was given)
angleRad = deg2rad(launchAngle);

% position arrays
x = zeros(1, n);% horizontal position
y = zeros(1, n); % vertical position

%velocity arrays  
vx = zeros(1, n);% horizontal velocity
vy = zeros(1, n);% vertical velocity

% Start at origin with zero velocity
x(1) = 0;
y(1) = 0;
vx(1) = 0; %horizontal component of velocity
vy(1) = 0; %vertical component of velocity

%% cross-sectional area
area = pi * (diameter/2)^2; % m^2: area of a circle equation

fprintf('Initial conditions:\n');
fprintf('Launch angle: %.1f degrees from horizontal\n', launchAngle);
fprintf('Thrust: %.0f N for %.1f seconds\n', thrust, burnTime);
fprintf('Mass: %.0f kg\n', mass);

%% simulation loop
for i = 2:n %for every data point in the array
    % Current time
    time = t(i);
    % Current speed (velocity mag)
    speed = sqrt(vx(i-1)^2 + vy(i-1)^2); %equation for speed
    
    %% Forces
    
    % thrust (only during burn time)
    if time <= burnTime % only burns for 5 sec
        thrustX = thrust * cos(angleRad);   % horizontal component
        thrustY = thrust * sin(angleRad);   % vertical component
    else
        thrustX = 0;
        thrustY = 0;
    end
    
    % drag(opposes motion, always present if velocity is nonzero)
    dragMag = 0.5 * rho * speed^2 * area * dragCoeff; %given equation
    if speed > 0
        %opposes velocity direction
        dragX = -dragMag * (vx(i-1) / speed); %dividing the previous velocity data pt by speed normalizes the vector
        dragY = -dragMag * (vy(i-1) / speed);
    else
        dragX = 0;
        dragY = 0;
    end
    
    % gravity (always present)
    gravityX = 0; %no gravity in the x direction
    gravityY = -mass * g; % negative bc downward
    
    % Net force 
    totalFx = thrustX + dragX + gravityX; %negative signs already accounted for 
    totalFy = thrustY + dragY + gravityY;
    
    %% accelerations (F = ma)
    ax = totalFx / mass;
    ay = totalFy / mass;
    
    %% update velocities (v = v0 + a*dt)
    vx(i) = vx(i-1) + ax * dt; 
    vy(i) = vy(i-1) + ay * dt;
    
    %% update positions using Euler method (numerical integration)
    x(i) = x(i-1) + vx(i-1)*dt + 0.5*ax*dt^2; %old position + displacement after 0.1 sec (1st order correction from velocity) + displacement after 0.1 sec (2nd order correction from acceleration)
    y(i) = y(i-1) + vy(i-1)*dt + 0.5*ay*dt^2; %same as above for y
    
    %% stop if rocket hits ground
    if y(i) <= 0 && i > 10  % small delay to avoid issues at beginning when y(i) = 0
        y(i) = 0; % set at ground level
        % cut arrays to actual flight time rather than what i put at beginning
        t = t(1:i);
        x = x(1:i);
        y = y(1:i);
        vx = vx(1:i);
        vy = vy(1:i);
        break;
    end
end

%% stats
[maxHeight, maxIdx] = max(y); % apogee is max height, maxIdx will give the time
apogeeTime = t(maxIdx);
apogeeX = x(maxIdx); %x coord at apogee
totalRange = x(end); %how far rocket went horizontally
flightTime = t(end);%total flight time
speed = sqrt(vx.^2 + vy.^2);% total speed array

%% show results
fprintf('\n Simulation Results \n');
fprintf('Launch Angle: %.1f degrees from horizontal\n', launchAngle); % %.1f means number will be displayed w 1 digit after decimal point
fprintf('Apogee: %.1f m at t=%.1f s\n', maxHeight, apogeeTime); % %d for integers, %f for floating-point numbers, %s for strings
fprintf('Range: %.1f m\n', totalRange);
fprintf('Flight Time: %.1f s\n', flightTime);
fprintf('Max Speed: %.1f m/s\n', max(speed));

%% close existing figures and make new ones
close all; 

%% Trajectory (x v y)
figure(1);
plot(x, y, 'b-', 'LineWidth', 2); % b- is blue line
hold on; %allows to add more plots to same figure without erasing existing plot (keeps plot active)
plot(apogeeX, maxHeight, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); %ro is red circle, 
plot(0, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); %go is green circle at origin (0,0)
xlabel('Horizontal Distance (m)');
ylabel('Altitude (m)');
title('Rocket Trajectory');
legend('Flight Path', 'Apogee', 'Launch', 'Location', 'best'); %best puts legend in optimal location to not block data
grid on;
axis equal; %keeps proportions of x and y axis equal

%makes matlab draw the figure
drawnow;

%% Velocity vs time  
figure(2);
subplot(3,1,1); %creates first subplot in 3 row by 1 column layout (this is first graph)
plot(t, vx, 'r-', 'LineWidth', 1.5); %horizontal v versus time red line
xlabel('Time (s)');
ylabel('Horizontal Velocity (m/s)');
title('Velocity Components vs Time');
grid on;

subplot(3,1,2);
plot(t, vy, 'b-', 'LineWidth', 1.5); %vertical velocity versus time blue line
xlabel('Time (s)');
ylabel('Vertical Velocity (m/s)');
grid on;

subplot(3,1,3);
plot(t, speed, 'k-', 'LineWidth', 1.5); %speed v time black line
xlabel('Time (s)');
ylabel('Total Speed (m/s)');
grid on;

drawnow;

%% force vs time
figure(3);

% recalculating forces for plot (same as earlier)
thrustForce = zeros(size(t));
dragForce = zeros(size(t));
gravityForce = mass * g * ones(size(t));

for i = 1:length(t)
    if t(i) <= burnTime
        thrustForce(i) = thrust;
    end
    currentSpeed = sqrt(vx(i)^2 + vy(i)^2);
    dragForce(i) = 0.5 * rho * currentSpeed^2 * area * dragCoeff;
end

plot(t, thrustForce, 'r-', 'LineWidth', 2); %thrust v time red line
hold on;
plot(t, dragForce, 'b-', 'LineWidth', 2); %drag v time blue line
plot(t, gravityForce, 'g-', 'LineWidth', 2); %gravity v time green line (constant)
xlabel('Time (s)'); 
ylabel('Force (N)');
title('Forces vs Time');
legend('Thrust', 'Drag', 'Weight', 'Location', 'best');
grid on;

drawnow;

%% Rocket animation using given code

fprintf('\nGenerating rocket animation...\n');
% For the animation, we need rotation angles
% Calculate actual rocket angle based on velocity direction
rocketAngle = zeros(size(t));
for i = 1:length(t) %for every time step
    if speed(i) > 0.1  % avoids division by zero (becuase first i will be 0)
        rocketAngle(i) = atan2(vy(i), vx(i)); %computes angle based on vertical and horizontal velocity components
    else
        rocketAngle(i) = angleRad; % use launch angle if no significant velocity
    end
end

% position array for animation function
positionArray = [x', y'];  % make into column vectors to plot

% rotation function (given)
rotating(positionArray, rocketAngle);

fprintf('Animation complete!\n');

%% given code for rotation
function rotating(pos, rot)
% Description:
% Animates rocket flight using provided algorithm
% Inputs:
% pos - position matrix with x and y [m]  
% rot - rotation angle array [rad]

    % Initialize figure for animation
    figure(4) %changed to figure 4 bc i already had figures 1-3
    clf; 
    
    x = pos(:,1);  % x positions [m]
    y = pos(:,2);  % y positions [m]
    
    w = 20;  % rocket width [m] (made larger because I couldnt see the rocket in the simulation)
    h = 60;  % rocket height [m]
    
    pgon = polyshape([-w/2,-w/2,0,w/2,w/2],[0,h,1.2*h,h,0]);  % rocket shape
    
    % Loop through position and rotation arrays to animate the rocket
    for i = 1:5:length(pos)  % Skip some frames for speed
        clf;  % Clear figure
        
        % Plot trajectory up to current point FIRST
        hold on;
        plot(x(1:i), y(1:i), 'r--', 'LineWidth', 1);
        
        % Update rocket position and rotation here AFTER trajectory
        pgonDraw = rotate(pgon, rad2deg(rot(i)-pi/2));  % rotate rocket [deg]
        pgonDraw = translate(pgonDraw, [x(i),y(i)]);    % move to position [m]
        plot(pgonDraw, 'FaceColor', 'blue');
        
        axis equal
        xlabel('Downrange [m]')
        ylabel('Altitude [m]')
        title(sprintf('Rocket Animation - Time: %.1f s', (i-1)*0.01*5))  % corrected time display
        grid on;
        
        % Set reasonable axis limits
        xlim([min(x)-100, max(x)+100]);
        ylim([-50, max(y)+100]);
    
        drawnow;  % Force MATLAB to update display
        pause(0.01)  % animation delay [s]
    
    end
end