%% Combined Method: Variation of Mechanical Trail & Scrub Radius Using a 3‐Point Wheel Plane
clc; clear; close all;

%% 1. Define Ball Joint Coordinates and Steering Axis
% Lower and upper ball joints (in meters)
P_upper = [-0.0134, -0.542, 0.357];    % Upper ball joint
P_lower = [-0.0064, -0.556, 0.166];   % Lower ball joint

% Compute the steering (kingpin) axis unit vector
axis_vec = P_upper - P_lower;
axis_unit = axis_vec / norm(axis_vec);

% Compute ground intersection (z = 0) of the steering axis
t_ground = -P_lower(3) / axis_unit(3);
P_k = P_lower + t_ground * axis_unit;  % This point (in x,y) serves as the reference

%% 2. Define the Initial Wheel Plane via Three Points
% These three points define the wheel plane in its neutral configuration.
% Replace these example points with your actual measurements.
P1 = [0, 1, 0];     
P2 = [0.5, 1, 0];   
P3 = [0, 1, 0.5];   

% Compute the initial plane normal (n0)
v1 = P2 - P1;
v2 = P3 - P1;
n0 = cross(v1, v2);
n0 = n0 / norm(n0);

%% 3. Define Hub Center and Tire Parameters
R_tire = 0.228;  % Tire radius (m)
% For this example, choose the hub center as the horizontal average of the ball joints 
% and set its height to the tire radius.
P_h = [(P_lower(1)+P_upper(1))/2, (P_lower(2)+P_upper(2))/2, R_tire];

%% 4. Set Up Steering Angle Variation
theta_max = deg2rad(30);  % Maximum steering angle (e.g., 30°) in radians
nAngles = 100;
steerAngles = linspace(-theta_max, theta_max, nAngles);

% Preallocate arrays for mechanical trail and scrub radius
mechanical_trail = zeros(size(steerAngles));
scrub_radius = zeros(size(steerAngles));

%% 5. Loop Over Steering Angles
for i = 1:length(steerAngles)
    angle = steerAngles(i);
    
    % Compute rotation matrix about the steering axis using Rodrigues’ formula
    R = rotationMatrix(axis_unit, angle);
    
    % Rotate the hub center about the steering axis (using P_lower as pivot)
    P_h_rot = P_lower + (R * (P_h - P_lower)')';
    
    % Rotate the initial wheel plane normal to obtain the new wheel plane normal
    n_new = (R * n0')';
    
    % To define the tire circle in the rotated wheel plane, determine two orthonormal 
    % vectors (u_in_plane and v_in_plane) spanning the plane.
    temp = [1, 0, 0];
    u_in_plane = temp - dot(temp, n_new)*n_new;
    if norm(u_in_plane) < 1e-6
        temp = [0, 1, 0];
        u_in_plane = temp - dot(temp, n_new)*n_new;
    end
    u_in_plane = u_in_plane / norm(u_in_plane);
    v_in_plane = cross(n_new, u_in_plane);
    v_in_plane = v_in_plane / norm(v_in_plane);
    
    % Parameterize the tire circle:
    % P(phi) = P_h_rot + R_tire*(cos(phi)*u_in_plane + sin(phi)*v_in_plane)
    % The lowest (contact) point is the one with minimum z. Because the circle is
    % smooth, the minimum is achieved when [cos(phi); sin(phi)] is opposite to the 
    % vector formed by the z‑components of u_in_plane and v_in_plane.
    %
    % Compute phi_min such that:
    % [cos(phi_min); sin(phi_min)] = -[u_in_plane(3); v_in_plane(3)]/norm([u_in_plane(3); v_in_plane(3)])
    phi_min = atan2(v_in_plane(3), u_in_plane(3)) + pi;
    
    % Compute the contact patch position (global coordinates)
    P_contact = P_h_rot + R_tire*(cos(phi_min)*u_in_plane + sin(phi_min)*v_in_plane);
    
    % Compute mechanical trail as the longitudinal (x-direction) offset between
    % the contact patch and the steering axis ground intersection.
    mechanical_trail(i) = P_contact(1) - P_k(1);
    
    % Compute scrub radius as the lateral (y-direction) offset.
    scrub_radius(i) = P_contact(2) - P_k(2);
end

%% 6. Plot the Results
figure;
subplot(2,1,1);
plot(rad2deg(steerAngles), mechanical_trail, 'b-', 'LineWidth', 2);
xlabel('Steering Angle (deg)');
ylabel('Mechanical Trail (m)');
title('Variation of Mechanical Trail with Steering Angle');
grid on;

subplot(2,1,2);
plot(rad2deg(steerAngles), scrub_radius, 'r-', 'LineWidth', 2);
xlabel('Steering Angle (deg)');
ylabel('Scrub Radius (m)');
title('Variation of Scrub Radius with Steering Angle');
grid on;

%% Helper Function: Rotation Matrix About an Arbitrary Axis
function R = rotationMatrix(u, theta)
    % u must be a unit vector (axis of rotation)
    % theta is the rotation angle in radians.
    %
    % Rodrigues' rotation formula:
    %   R = I + sin(theta)*K + (1-cos(theta))*(K^2)
    % where K is the skew-symmetric matrix of u.
    
    K = [   0    -u(3)  u(2);
          u(3)    0    -u(1);
         -u(2)  u(1)    0];
    R = eye(3) + sin(theta)*K + (1-cos(theta))*(K*K);
end
