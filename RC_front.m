%% FSAE Car Roll Centre and Jacking Force Simulation
clear; clc; close all;

%% Constants and Initial Values
W = 290*9.8;          % Weight of the vehicle (N)
rg = 1.5;             % Roll gradient (unused in current code)
tR = 1.16;            % Rear track width (m)
H = 0.259;            % Height of the CG (m)
a = 1.53/2;           % Distance from CG to front axle (m)
l = 1.53;             % Wheelbase (m)
motion_ratio_rear = 2.09;
Mx = 0;               % Moment from tire forces (currently zero; update if needed)
wheelcentre_L = 0.5;  % Lateral location for left wheel (assumed contact point)
wheelcentre_R = 0.5;  % Lateral location for right wheel (assumed contact point)
M_nsm = 50;           % Mass used for NSM inertial force (kg)
K = 108302.6714;      % Rear wheel spring rate (N/m)

% Suspension Geometry Initialization
% Define inboard (fixed, chassis-mounted) points for rear wheels (these remain fixed now)
inboard_RR_init = [0.255, 0.29502; 0.205, 0.13809];
inboard_LR_init = [-0.255, 0.29502; -0.205, 0.13809];

% Define outboard points for rear wheels (these will be rotated)
outboard_RR_init = [0.543814, 0.35235; 0.55603, 0.16335];
outboard_LR_init = [-0.543814, 0.35235; -0.55603, 0.16335];

% Rear wheelbase contact points at ground level (initial)
rear_base_RR_init = [0.580, 0];  
rear_base_LR_init = [-0.580, 0];

%% Define Range of Roll Angles
roll_values_deg = 0:0.5:20;  
n_roll = length(roll_values_deg);

% Preallocate arrays to store results
x_roll_centres = zeros(n_roll, 1);
y_roll_centres = zeros(n_roll, 1);
JF_R = zeros(n_roll, 1);
JF_L = zeros(n_roll, 1);

%% Main Loop: Simulate for Each Roll Angle
for i = 1:n_roll
    % Current roll angle in degrees and radians
    roll_deg = roll_values_deg(i);
    roll_rad = deg2rad(roll_deg);
    
    % Reinitialize geometry for this iteration to avoid cumulative updates
    rear_base_RR = rear_base_RR_init;
    rear_base_LR = rear_base_LR_init;
    
    % Fix the inboard points (vehicle/chassis remains fixed)
    inboard_RR = inboard_RR_init;
    inboard_LR = inboard_LR_init;
    
    % Rotate the outboard points instead
    Rmat = [cos(roll_rad) -sin(roll_rad); sin(roll_rad) cos(roll_rad)];
    outboard_RR = (Rmat * outboard_RR_init')';
    outboard_LR = (Rmat * outboard_LR_init')';
    
    %% Calculate Instantaneous Centers (IC) for Each Rear Wheel
    % Rear Left Wheel
    m_LR1 = (outboard_LR(1,2) - inboard_LR(1,2)) / (outboard_LR(1,1) - inboard_LR(1,1));
    b_LR1 = inboard_LR(1,2) - m_LR1 * inboard_LR(1,1);
    
    m_LR2 = (outboard_LR(2,2) - inboard_LR(2,2)) / (outboard_LR(2,1) - inboard_LR(2,1));
    b_LR2 = inboard_LR(2,2) - m_LR2 * inboard_LR(2,1);
    
    x_IC1 = (b_LR2 - b_LR1) / (m_LR1 - m_LR2);
    y_IC1 = m_LR1 * x_IC1 + b_LR1;
    
    % Rear Right Wheel
    m_RR1 = (outboard_RR(1,2) - inboard_RR(1,2)) / (outboard_RR(1,1) - inboard_RR(1,1));
    b_RR1 = inboard_RR(1,2) - m_RR1 * inboard_RR(1,1);
    
    m_RR2 = (outboard_RR(2,2) - inboard_RR(2,2)) / (outboard_RR(2,1) - inboard_RR(2,1));
    b_RR2 = inboard_RR(2,2) - m_RR2 * inboard_RR(2,1);
    
    x_IC2 = (b_RR2 - b_RR1) / (m_RR1 - m_RR2);
    y_IC2 = m_RR1 * x_IC2 + b_RR1;
    
    %% Calculate the Roll Center by Intersecting Lines from Each IC to Their Respective Wheel Contact Points
    m_LR_IC = (rear_base_LR(2) - y_IC1) / (rear_base_LR(1) - x_IC1);
    b_LR_IC = y_IC1 - m_LR_IC * x_IC1;
    
    m_RR_IC = (rear_base_RR(2) - y_IC2) / (rear_base_RR(1) - x_IC2);
    b_RR_IC = y_IC2 - m_RR_IC * x_IC2;
    
    x_roll_center = (b_RR_IC - b_LR_IC) / (m_LR_IC - m_RR_IC);
    y_roll_center = m_LR_IC * x_roll_center + b_LR_IC;
    
    % Store roll center positions
    x_roll_centres(i) = x_roll_center;
    y_roll_centres(i) = y_roll_center;
    
    %% Lateral Load Transfer Calculation
    delta_WR_Ay = W / tR * (H/2 + (a/l * y_roll_center));
    
    %% Jacking Force Calculation
    % Placeholder lateral forces from a tire model (update if nonzero values are available)
    Fy_L = 0;  % Lateral force on left tire (N)
    Fy_R = 0;  % Lateral force on right tire (N)
    
    theta_R = atan2(y_IC2, rear_base_RR(1) - x_IC2);
    theta_L = atan2(y_IC1, rear_base_LR(1) - x_IC1);
    
    jacking_force_L_LF = Fy_L * tan(theta_L);
    jacking_force_R_LF = -Fy_R * tan(theta_R);
    
    jacking_force_R_Mx = Mx / (rear_base_RR(1) - x_IC2);
    jacking_force_L_Mx = Mx / (rear_base_LR(1) - x_IC1);
    
    alpha_R = atan2(wheelcentre_R - y_IC2, rear_base_RR(1) - x_IC2);
    alpha_L = atan2(wheelcentre_L - y_IC1, rear_base_LR(1) - x_IC1);
    jacking_force_R_NSM = -M_nsm * 9.81 * a * alpha_R;
    jacking_force_L_NSM = M_nsm * 9.81 * a * alpha_L;
    
    JF_L(i) = jacking_force_L_NSM + jacking_force_L_Mx + jacking_force_L_LF;
    JF_R(i) = jacking_force_R_NSM + jacking_force_R_Mx + jacking_force_R_LF;
    
    %% Vertical Displacement due to Load Transfer
    spring_displacement = delta_WR_Ay / K;
    vertical_displacement = motion_ratio_rear * spring_displacement;
    
    rear_base_LR(2) = rear_base_LR_init(2) + vertical_displacement;
    rear_base_RR(2) = rear_base_RR_init(2) - vertical_displacement;
end

%% Plotting Results
figure;
plot(roll_values_deg, y_roll_centres*1000, 'b-o');
xlabel('Roll Angle (deg)');
ylabel('Roll Center Y-Coordinate (mm)');
title('Variation of Roll Center Y-Coordinate with Roll Angle');
grid on;

figure;
subplot(2,1,1);
plot(roll_values_deg, JF_L, 'r-o');
xlabel('Roll Angle (deg)');
ylabel('Jacking Force (N)');
title('Jacking Force on Left Wheel vs. Roll Angle');
grid on;

subplot(2,1,2);
plot(roll_values_deg, JF_R, 'b-o');
xlabel('Roll Angle (deg)');
ylabel('Jacking Force (N)');
title('Jacking Force on Right Wheel vs. Roll Angle');
grid on;
