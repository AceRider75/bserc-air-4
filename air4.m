% Clear any old variables and close plots
clear all;
close all;
clc;


Cr = 0.9;          % Root chord (m)
Ct = 0.15;         % Tip chord (m)
b = 1.5;           % Span (m)
e = 0.89;          % Oswald efficiency factor
CD0 = 0.03;        % Zero-lift drag coeff
Cm0 = 0.01;        % Zero-alpha pitching moment coeff
CL_alpha = 2.92;   % Lift slope (per rad)
CL_deltae = 0.265; % Lift due to elevator (per rad? assuming, but document doesn't specify units—check if needed)
Cm_alpha = -0.292; % Pitching moment slope (per rad)
Cm_deltae = -0.4;  % Pitching moment due to elevator (per rad? same note)
m = 3.5;           % Mass (kg)
rho = 1.225;       % Air density (kg/m^3)
g = 10;            % Gravity (m/s^2)
CL0 = 0;           % Assume zero-lift CL (not given, standard assumption)

W = m * g;         % Weight (N)


lambda = Ct / Cr;                    % Taper ratio
S = (b / 2) * Cr * (1 + lambda);    % Planform area (m^2) — note: document has S = b/2 * Cr (1+λ), but standard is (Cr+Ct)*b/2 = Cr*(1+λ)*b/2, yes same
AR = b^2 / S;                        % Aspect ratio
k = 1 / (pi * e * AR);              % Induced drag factor

% Display preprocessing results 
disp('Preprocessing Results:');
disp(['Taper Ratio (lambda): ', num2str(lambda)]);
disp(['Planform Area (S): ', num2str(S), ' m^2']);
disp(['Aspect Ratio (AR): ', num2str(AR)]);
disp(['Induced Drag Factor (k): ', num2str(k)]);


alpha_deg = 0:0.5:12;                % From 0 to 12 deg, step 0.5
num_points = length(alpha_deg);    


alpha_rad = zeros(1, num_points);
deltae_trim = zeros(1, num_points);
CL_trim = zeros(1, num_points);
V = zeros(1, num_points);
CD = zeros(1, num_points);
Tr = zeros(1, num_points);
Pr = zeros(1, num_points);
L_over_D = zeros(1, num_points);
L32_over_D = zeros(1, num_points);


for i = 1:num_points
    alpha_rad(i) = deg2rad(alpha_deg(i));  % Convert to radians
    
    % Trim elevator deflection from Cm = 0
    deltae_trim(i) = - (Cm0 + Cm_alpha * alpha_rad(i)) / Cm_deltae;
    
    % Corresponding lift coefficient
    CL_trim(i) = CL0 + CL_alpha * alpha_rad(i) + CL_deltae * deltae_trim(i);
    
    % Flight velocity from lift = weight
    V(i) = sqrt(2 * W / (rho * S * CL_trim(i)));
    
    % Drag coefficient
    CD(i) = CD0 + k * CL_trim(i)^2;
    
    % Thrust required
    Tr(i) = 0.5 * rho * V(i)^2 * S * CD(i);
    
    % Power required
    Pr(i) = Tr(i) * V(i);
    
    % Performance metrics
    L_over_D(i) = CL_trim(i) / CD(i);
    L32_over_D(i) = CL_trim(i)^(3/2) / CD(i);
end


% Plot vs angle of attack (alpha_deg)
figure(1); 
subplot(4,2,1); plot(alpha_deg, deltae_trim); xlabel('Alpha (deg)'); ylabel('Delta_e (rad)'); title('Trim Elevator Deflection');
subplot(4,2,2); plot(alpha_deg, CL_trim); xlabel('Alpha (deg)'); ylabel('CL'); title('Lift Coefficient');
subplot(4,2,3); plot(alpha_deg, CD); xlabel('Alpha (deg)'); ylabel('CD'); title('Drag Coefficient');
subplot(4,2,4); plot(alpha_deg, Tr); xlabel('Alpha (deg)'); ylabel('Thrust (N)'); title('Thrust Required');
subplot(4,2,5); plot(alpha_deg, Pr); xlabel('Alpha (deg)'); ylabel('Power (W)'); title('Power Required');
subplot(4,2,6); plot(alpha_deg, V); xlabel('Alpha (deg)'); ylabel('V (m/s)'); title('Flight Velocity');
subplot(4,2,7); plot(alpha_deg, L_over_D); xlabel('Alpha (deg)'); ylabel('CL/CD'); title('CL/CD');
subplot(4,2,8); plot(alpha_deg, L32_over_D); xlabel('Alpha (deg)'); ylabel('CL^{3/2}/CD'); title('CL^{3/2}/CD');

% Also plot vs velocity (V) as alternative
figure(2);
subplot(4,2,1); plot(V, deltae_trim); xlabel('V (m/s)'); ylabel('Delta_e (rad)'); title('Trim Elevator Deflection');
subplot(4,2,2); plot(V, CL_trim); xlabel('V (m/s)'); ylabel('CL'); title('Lift Coefficient');
subplot(4,2,3); plot(V, CD); xlabel('V (m/s)'); ylabel('CD'); title('Drag Coefficient');
subplot(4,2,4); plot(V, Tr); xlabel('V (m/s)'); ylabel('Thrust (N)'); title('Thrust Required');
subplot(4,2,5); plot(V, Pr); xlabel('V (m/s)'); ylabel('Power (W)'); title('Power Required');
subplot(4,2,6); plot(V, alpha_deg); xlabel('V (m/s)'); ylabel('Alpha (deg)'); title('Trim Angle of Attack');
subplot(4,2,7); plot(V, L_over_D); xlabel('V (m/s)'); ylabel('CL/CD'); title('CL/CD');
subplot(4,2,8); plot(V, L32_over_D); xlabel('V (m/s)'); ylabel('CL^{3/2}/CD'); title('CL^{3/2}/CD');


results_table = table(alpha_deg', V', deltae_trim', CL_trim', CD', Tr', Pr', L_over_D', L32_over_D', ...
    'VariableNames', {'Alpha_deg', 'V_mps', 'Deltae_rad', 'CL', 'CD', 'Thrust_N', 'Power_W', 'CL_CD', 'CL32_CD'});
disp(results_table);  % Show in Command Window
writetable(results_table, 'UAV_Results.xlsx');
