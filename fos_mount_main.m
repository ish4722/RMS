%--- Tire FOS Analysis with Load Transfer for 6 Rods and Three Cases ---
clear; clc;

%--- 0) Initialize Tire Model ---
tireDataFilePath = 'Round_6_Hoosier_R25B_20p5x7_13_on_7in_10psi_PAC02_UM4.tir';
tireModel = TireModel_MINE(tireDataFilePath);
Vx = 11.1760;      % vehicle forward speed [m/s]

%--- 1) Vehicle & CG parameters ---
m    = 290;          % mass [kg]
g    = 9.81;        % gravity [m/s^2]
h_cg = 0.259;        % CG height [m]
L    = 1.53;        % wheelbase [m]
l_f  = 1.53/2;        % CG to front axle [m]
l_r  = 1.53/2;    % CG to rear axle [m]
t    = 1.16;        % track width [m]

%--- 2) Static wheel loads ---
Fz_static_front = (m*g)*(l_r/L);
Fz_static_rear  = (m*g)*(l_f/L);
Fz0_vec = [Fz_static_front/2, Fz_static_front/2, Fz_static_rear/2, Fz_static_rear/2];

names = {'fl','fr','rl','rr'};
peakB = struct(); peakC = struct(); peakA = struct();

%--- (A) Braking Case ---
alpha0 = 0;
kappaRange = linspace(-0.25, 0.25, 50);
Fx_nom = zeros(1,4);
Mz_nom = zeros(1,4);
My_nom = zeros(1,4);
for k = 1:4
    % per-wheel longitudinal responses over kappa
    Fx_vals = tireModel.longc(Fz0_vec(k), alpha0, 0, kappaRange);
    Mz_vals = tireModel.Alignc(Fz0_vec(k), alpha0, 0, kappaRange);
    My_vals = tireModel.Rolling(Fz0_vec(k), 0, kappaRange, Vx);
    Fx_nom(k) = max(abs(Fx_vals));
    Mz_nom(k) = max(abs(Mz_vals));
    My_nom(k) = max(abs(My_vals));
end
% total peak braking force and deceleration
totFxPeak = sum(Fx_nom);
a_x = totFxPeak / m;
dFz_long = m * a_x * h_cg / L;
FzB = [Fz0_vec(1)+dFz_long/2, Fz0_vec(2)+dFz_long/2, Fz0_vec(3)-dFz_long/2, Fz0_vec(4)-dFz_long/2];
for k = 1:4
    peakB.(names{k}) = [Fx_nom(k), 0, Mz_nom(k), 0, My_nom(k)];
end

%--- (B) Cornering Case ---
alphaRange = linspace(-0.2618, 0.2618, 50);
kappa0 = 0;
Fy_nom = zeros(1,4);
Fx_c = zeros(1,4);
Mz_c = zeros(1,4);
Mx_c = zeros(1,4);
for k = 1:4
    Fy_vals = tireModel.latc(Fz0_vec(k), 0, 0, 0);
    Fx_vals = zeros(size(alphaRange));
    Mz_vals = zeros(size(alphaRange));
    Mx_vals = zeros(size(alphaRange));
    Fy_vals = zeros(size(alphaRange));
    for i = 1:length(alphaRange)
        al = alphaRange(i);
        Fy_vals(i) = tireModel.latc(Fz0_vec(k), al, 0, kappa0);
        Fx_vals(i) = tireModel.longc(Fz0_vec(k), al, 0, kappa0);
        Mz_vals(i) = tireModel.Alignc(Fz0_vec(k), al, 0, kappa0);
        Mx_vals(i) = tireModel.Overturn(Fz0_vec(k), al, 0, kappa0);
    end
    Fy_nom(k) = max(abs(Fy_vals));
    Fx_c(k)  = max(abs(Fx_vals));
    Mz_c(k)  = max(abs(Mz_vals));
    Mx_c(k)  = max(abs(Mx_vals));
end
% total peak lateral force and lateral acceleration
totFyPeak = sum(Fy_nom);
a_y = totFyPeak / m;
dFz_lat = m * a_y * h_cg / t;
FzC = [Fz0_vec(1)-dFz_lat/2, Fz0_vec(2)+dFz_lat/2, Fz0_vec(3)-dFz_lat/2, Fz0_vec(4)+dFz_lat/2];
for k = 1:4
    peakC.(names{k}) = [Fx_c(k), Fy_nom(k), Mz_c(k), Mx_c(k), 0];
end

%--- (C) Combined Case ---
for k = 1:4
    pB = peakB.(names{k});
    pC = peakC.(names{k});
    peakA.(names{k}) = max([pB; pC], [], 1);
end

%--- 4) Assemble results table ---
cases = {'Braking','Cornering','Combined'};
rows = numel(cases)*4*6;
Case = cell(rows,1); Wheel = cell(rows,1); Rod = zeros(rows,1);
BracketStress_MPa = zeros(rows,1); BracketFOS = zeros(rows,1);
WeldStress_MPa = zeros(rows,1); WeldFOS = zeros(rows,1);
row = 0;
for ci = 1:numel(cases)
    caseName = cases{ci};
    peakStruct = eval(['peak' caseName(1)]);
    for wi = 1:4
        w = names{wi};
        loads = peakStruct.(w);
        [~, sb, fb, sw, fw] = calculateFOSFull(loads(1), loads(2), loads(3), loads(4), loads(5), Fz0_vec(wi));
        for ri = 1:6
            row = row + 1;
            Case{row} = caseName;
            Wheel{row} = upper(w);
            Rod(row) = ri;
            BracketStress_MPa(row) = sb(ri);
            BracketFOS(row) = fb(ri);
            WeldStress_MPa(row) = sw(ri);
            WeldFOS(row) = fw(ri);
        end
    end
end
resultsTable = table(Case, Wheel, Rod, BracketStress_MPa, BracketFOS, WeldStress_MPa, WeldFOS);

%--- 5) Write to Excel ---
writetable(resultsTable, 'TireFOS_AllCases.xlsx');
fprintf('All cases FOS saved to TireFOS_AllCases.xlsx\n');

%% calculateFOSFull: returns rod stresses & FOS for 6 mounts
function [Frod, stressB, fosB, stressW, fosW] = calculateFOSFull(Fx, Fy, Mz, Mx, My, Fz)
    syms F1 F2 F3 F4 F5 F6;
    eq1 = F1*0.40 + F2*(-0.59) + F3*0.49 + F4*(-0.56) + F5*0.99 + 0.69*F6 == Fx;
    eq2 = -(F1*0.09) - (F2*0.02) - (F3*0.16) - (F4*0.18) + F5*(-0.11) + 0.72*F6 == Fy;
    eq3 = F1*(-0.91) - F2*0.80 + F3*(-0.85) - F4*0.80 == Fz;
    eq4 = F1*(-0.91)*0.124 + F2*(-0.80)*0.124 + F3*(-0.85)*0.244 + F4*(-0.80)*0.244 == Mx;
    eq5 = F1*0.40*0.807 + F2*(-0.59)*0.807 + F3*(-0.16)*0.788 + F4*(-0.22)*0.788 + F5*0.99*0.15 + F6*0.69*0 == My;
    eq6 = F1*0.29*0.577 + F2*(-0.33)*0.577 + F3*0.36*0.564 + F4*(-0.33)*0.564 + F5*(-0.11)*0.975 + F6*(-0.72)*0.099 == Mz;
    sol = solve([eq1, eq2, eq3, eq4, eq5, eq6], [F1, F2, F3, F4, F5, F6]);
    F = double([sol.F1; sol.F2; sol.F3; sol.F4; sol.F5; sol.F6]);

    thickness = 2; hole_dia = 10; 
    width = 20; weld_area = thickness * width;
    bracket_area = 2 * (width-hole_dia);
    kt = 2.18; Sy_b = 305; Sy_w = 180;

    Frod = abs(F) / 2;
    stressB = (Frod ./ bracket_area) * kt;
    stressW = Frod ./ weld_area;
    fosB = Sy_b ./ stressB;
    fosW = Sy_w ./ stressW;
end
