% Given Forces and Moments
Fx = 3.572183053729892e+03; % N
Fy = 3.087041385355991e+03; %
Fz = 2.534422566371682e+03; % N
Mx = -66.461539154186870; % Nm
My = 0; % Nm
Mz = -1.965138290235688e+02; % Nm

% Direction vectors from the image
d1 = [711.60-894.37, 325.24-384.86, 970.10-1233.76];
d2 = [1062.39-894.37, 333.18-384.86, 944.26-1233.76];
d3 = [676.70-900.98, 188.21-195.86, 944.36-1245.98];
d4 = [1065.66-900.98, 158.72-195.86, 871.62-1245.98];
d5 = [1,5, 3]; % direction vector for push rod
d6 = [188.62, -188.62, 0]; % direction vector for tie rod

% Normalize direction vectors
d1 = d1/norm(d1);
d2 = d2/norm(d2);
d3 = d3/norm(d3);
d4 = d4/norm(d4);
d5 = d5/norm(d5);
d6 = d6/norm(d6);

% Position vectors from the center of the wheel to each rod's point of action
r1 = [711.60, 325.24, 970.10];
r2 = [1062.39, 333.18, 944.26];
r3 = [676.70, 188.21, 944.36];
r4 = [1065.66, 158.72, 871.62];
r5 = [500, 150, 900]; % position for push rod
r6 = [188.62, 0, 0]; % position for tie rod

% Calculate cross products
m1 = cross(r1, d1);
m2 = cross(r2, d2);
m3 = cross(r3, d3);
m4 = cross(r4, d4);
m5 = cross(r5, d5);
m6 = cross(r6, d6);

% Normalize moment arms
m1 = m1/norm(m1);
m2 = m2/norm(m2);
m3 = m3/norm(m3);
m4 = m4/norm(m4);
m5 = m5/norm(m5);
m6 = m6/norm(m6);

% Construct the matrix
A = [d1(1), d2(1), d3(1), d4(1), d5(1), d6(1);
     d1(2), d2(2), d3(2), d4(2), d5(2), d6(2);
     d1(3), d2(3), d3(3), d4(3), d5(3), d6(3);
     m1(1), m2(1), m3(1), m4(1), m5(1), m6(1);
     m1(2), m2(2), m3(2), m4(2), m5(2), m6(2);
     m1(3), m2(3), m3(3), m4(3), m5(3), m6(3)];

% Given Forces and Moments vector
B = [Fx; Fy; Fz; Mx; My; Mz];

% Solve for the forces in each rod
F = inv(A) * B;

% Display the results
fprintf('Forces in the rods:\n');
for i = 1:length(F)
    fprintf('F%d: %.2f N\n', i, F(i));
end
