%masses are taken from 50 to 125 at an interval of 25
%slip angles are taken from 0 to 2 at interval of 0.05

% Define constants
QBZ1 = 9.53877;
QBZ2 = -5.13125;
QBZ3 = -3.05044;
QBZ4 = 2.3916;
QBZ5 = -3.338;
QBZ9 = 0;
QBZ10 = 4.19999e-07;
QCZ1 = 1.04832;
QDZ1 = 0.142958;
QDZ2 = -0.0487339;
QDZ3 = 1.1459;
QDZ4 = -22.3675;
QDZ6 = -0.0103289;
QDZ7 = 0.0289023;
QDZ8 = -0.856;
QDZ9 = -0.1181;
QEZ1 = -4.20626;
QEZ2 = -4.3916;
QEZ3 = 0.209767;
QEZ4 = 0.291079;
QEZ5 = -10;
QHZ1 = -0.00670191;
QHZ2 = 0.00830384;
QHZ3 = -0.05;
QHZ4 = -0.0696;
UNLOADED_RADIUS = 0.2603;
VERTICAL_STIFFNESS = 109980;
VERTICAL_DAMPING = 205.7297;
lambda_t = 1;
lambda_Ky = 1;
lambda_mu_y = 0.67;
lamGamy = 1;
Fz0 = 1054;
Lfzo = 1;

% Define load range
Fz_range =50*9.8:1*9.8:125*9.8;
Fz_alpha = numel(Fz_range);

% Initialize figure
figure;
hold on;
colors = ['r', 'g', 'b', 'k', 'm']; % Colors for different lines

% Loop over different loads
for k = 1:Fz_alpha
    R0 = UNLOADED_RADIUS;
    Fzo1 = Fz0 * Lfzo;
    dfz = (Fz_range(k) - Fzo1) / Fzo1;
    gamma_z = 0; % Assume no camber
    alpha_range = 0:0.05:0.2;
    num_alpha = numel(alpha_range);
    AlTarray = zeros(1, num_alpha);
    SHt = QHZ1 + QHZ2 * dfz + QHZ3 * gamma_z + QHZ4 * dfz * gamma_z;

    for j = 1:num_alpha

        % Calculate the scaled slip angles
        alpha_t = alpha_range(j) + SHt;

        % Calculate the coefficients
        Dt = Fz_range(k) * (QDZ1 + QDZ2 * dfz) * (1 + QDZ3 * gamma_z + QDZ4 * gamma_z^2) * (R0 / Fz0) * lambda_t;
        Ct = QCZ1;
        Bt = (QBZ1 + QBZ2 * dfz + QBZ3 * dfz^2) * (1 + QBZ4 * gamma_z + QBZ5 * abs(gamma_z)) * (lambda_Ky / lambda_mu_y);
        Et = (QEZ1 + QEZ2 * dfz + QEZ3 * dfz^2) * (1 + (QEZ4 + QEZ5 * gamma_z) * (2 / pi) * atan(Bt * Ct * alpha_t));

        % Calculate the pneumatic trail
        t = Dt * cos(Ct * atan(Bt * alpha_t - Et * (Bt * alpha_t - atan(Bt * alpha_t)))) * cos(alpha_range(j));

        % Calculate lateral force using the latp function
        lamFz0 = 1;
        lamCy = 1;
        lamMuy = 0.67;
        lamEy = 1;
        lamKy = 1;
        lamHy = 1;
        lamVy = 1;
        lamGamy = 1;
    
        pCy1 = 1.51911;
        pDy1 = -2.6751;
        pDy2 = 0.571186;
        pDy3 = 2;
        pEy1 = -0.139324;
        pEy2 = 0.215641;
        pEy3 = 1;
        pEy4 = 1;
        pKy1 = -53.4851;
        pKy2 = 1.40751;
        pKy3 = 2.21073;
        pHy1 = -0.00135288;
        pHy2 = -0.000348694;
        pHy3 = 0.176518;
        pVy1 = -0.027448;
        pVy2 = 0.0587239;
        pVy3 = 2.37307;
        pVy4 = -2.1928;
        fz0 = 1054;
        lamFz0 = 1;
        fz01 = fz0 * lamFz0;
        dFz = (Fz_range(k) - fz01) / fz01;
        Gamy = gamma_z * lamGamy;
        SHy = (pHy1 + pHy2 * dFz) * lamHy + pHy3 * Gamy;
        SVy =Fz_range(k) * ((pVy1 + pVy2 * dFz) * lamVy + (pVy3 + pVy4 * dFz) * Gamy) * lamMuy;
        Aly = alpha_t + SHy;
        Cy = pCy1 * lamCy;
        Muy = (pDy1 + pDy2 * dFz) * (1 - pDy3 * Gamy^2) * lamMuy;
        Dy = Muy * Fz_range(k);
        Ey = (pEy1 + pEy2 * dFz) * (1 - (pEy3 + pEy4 * Gamy) * sign(Aly)) * lamEy;
        Ky = pKy1 * fz0 * sin(2 * atan(Fz_range(k) / (pKy2 * fz0 * lamFz0))) * (1 - pKy3 * abs(Gamy)) * lamFz0 * lamKy;
        By = Ky / (Cy * Dy);
        Fy0 = latp(Fz_range(k), alpha_t, gamma_z);

       % Residual aligning torque calculation
        Dr = Fz_range(k) * (QDZ6 + QDZ7 * dfz) * (1 + QDZ8 * gamma_z + QDZ9 * gamma_z * dfz) * R0;
        Br = QBZ9*lambda_Ky/lambda_mu_y+QBZ10*By*Cy;
        SHr = QHZ1 + QHZ2 * dfz + QHZ3 * gamma_z + QHZ4 * dfz * gamma_z;
        alpha_r = alpha_range(j) + SHr;
        Mzr = Dr * cos(atan(Br * alpha_r))*cos(alpha_range(j));
        
        % Calculate the aligning torque
        Mz0 = t * Fy0 - Mzr;

        % Display the results
        AlTarray(j) = Mz0;
    end

    % Plot results
    plot(alpha_range, AlTarray, 'DisplayName', sprintf('Fz = %d N', Fz_range(k)), 'Color', colors(k));
end

% Customize the plot
xlabel('Slip Angle');
ylabel('Aligning Moment');
title('Aligning Moment vs. Slip Angle for Different Masses');
legend('show');
hold off;





% Define lateral force function
function Fy0 = latp(Fz, Al, gamma_z)
    lamFz0 = 1;
    lamCy = 1;
    lamMuy = 0.67;
    lamEy = 1;
    lamKy = 1;
    lamHy = 1;
    lamVy = 1;
    lamGamy = 1;

    pCy1 = 1.51911;
    pDy1 = -2.6751;
    pDy2 = 0.571186;
    pDy3 = 2;
    pEy1 = -0.139324;
    pEy2 = 0.215641;
    pEy3 = 1;
    pEy4 = 1;
    pKy1 = -53.4851;
    pKy2 = 1.40751;
    pKy3 = 2.21073;
    pHy1 = -0.00135288;
    pHy2 = -0.000348694;
    pHy3 = 0.176518;
    pVy1 = -0.027448;
    pVy2 = 0.0587239;
    pVy3 = 2.37307;
    pVy4 = -2.1928;
    Fz0 = 1054;
    lamFz0 = 1;
    Fz01 = Fz0 * lamFz0;
    dfz = (Fz - Fz01) / Fz01;
    Gamy = gamma_z * lamGamy;
    SHy = (pHy1 + pHy2 * dfz) * lamHy + pHy3 * Gamy;
    SVy = Fz * ((pVy1 + pVy2 * dfz) * lamVy + (pVy3 + pVy4 * dfz) * Gamy) * lamMuy;
    Aly = Al + SHy;
    Cy = pCy1 * lamCy;
    Muy = (pDy1 + pDy2 * dfz) * (1 - pDy3 * Gamy^2) * lamMuy;
    Dy = Muy * Fz;
    Ey = (pEy1 + pEy2 * dfz) * (1 - (pEy3 + pEy4 * Gamy) * sign(Aly)) * lamEy;
    Ky = pKy1 * Fz0 * sin(2 * atan(Fz / (pKy2 * Fz0 * lamFz0))) * (1 - pKy3 * abs(Gamy)) * lamFz0 * lamKy;
    By = Ky / (Cy * Dy);
    Fy0 = -(Dy * sin(Cy * atan(By * Aly - Ey * (By * Aly - atan(By * Aly)))) + SVy);
end