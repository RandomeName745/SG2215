clc; clear all; close all;

%% Constants
Minf = 2;
AOA = 1.5*pi/180;              % angle of attack
gamma = 1.4;

%% Wing properties
r = 210;    %[-] radius for wing curvature
c = 70;     %[-] chord of the wing

x = [7 17 27 37 47 57 67];
phi = asin((c/2 - x)/r);  %[rad] inclination of the wing at each pressure tap

% only measurement values for lower surface
% p: lower surface (plus)
% m: upper surface(minus)
theta_p15 = phi + AOA;
theta_m15 = phi - AOA;

phi0 = asin((c/2)/r);  %[rad] inclination of the wing at each pressure tap
theta0_p15 = phi0 + AOA;
theta0_m15 = phi0 - AOA;



%% experimental data
ptaps_sens = [22.22	22.311	22.27	22.304	22.234	22.464	22.518	22.63]; %[kPa/V] sensitivity
ptaps_offset = [0.202	0.202	0.204	0.204	0.199	0.197	0.203	0.202]; %[V] offset

ptaps_p15 = [48.129	41.635	35.315	30.077	25.953	22.554	19.322	41.244];    %[V] measured voltage
ptaps_m15 = [40.988	34.778	29.776	25.167	21.542	18.447	18.406	49.229]; %[V] measured voltage

p0_volt_p15 = 1.546;    %[V]
p0_volt_m15 = 1.556;    %[V]
K1 = 79.6737;   %[kPa/V] calibration constant
K0 = 0.0202;    %[kPa] colibratuib constant

%stagnation pressure 
patm = 101.2;   %[kPa] atmospheric pressure
p0_p15 = K1*p0_volt_p15 + K0 + patm;
p0_m15 = K1*p0_volt_m15 + K0 + patm;

%infinity pressure
pinf_p15 = p0_p15/((1+(gamma-1)/2*Minf^2)^(gamma/(gamma-1)));
pinf_m15 = p0_m15/((1+(gamma-1)/2*Minf^2)^(gamma/(gamma-1)));

%cp
cp_p15 = 2/(gamma*Minf^2)*(ptaps_p15./pinf_p15-1);
cp_m15 = 2/(gamma*Minf^2)*(ptaps_m15./pinf_m15-1);


%% linearized theory in comparison
cp_p15_lin = 2*theta_p15/(sqrt(Minf^2-1));
cp_m15_lin = 2*theta_m15/(sqrt(Minf^2-1));

%% exact shock wave theory (theta-beta-M-rel for oblique shock and Prandtl-Meyer for expansion)

for jj = 1:2
    
    if jj == 1
        theta = theta_p15;
        theta0 = theta0_p15;
%         pinf = pinf_p15;
    else
        theta = theta_m15;
        theta0 = theta0_m15;
%         pinf = pinf_m15;
    end
    
    beta_init = [0, pi/2];
    M_init = [0, 5];
    tol = 1e-6;

    
    % Calculate shock angle for oblique shock at leading edge
    beta0 = BisectionAlgorithm_ThetaBetaM(theta0, Minf, beta_init(1), beta_init(2), gamma, tol);
    % Calculate normal Mach number
    Mn1 = Minf*sin(beta0);
    % Calculate pressure ratio across oblqiue shock from normal shock relation
    p2_pinf = 1+(2*gamma)/(gamma+1)*(Mn1^2-1);
    % Calculate normal Mach number behind oblique shock
    Mn2 = sqrt((1 + (gamma - 1)/2 * Mn1^2) / (gamma*Mn1^2 - (gamma - 1)/2));

    % Calculate Mach number behind oblique shock
    M2 = Mn2 / sin(beta0 - theta0);
    

    nu0 = CalculatePrandtlMayerAngle(M2, gamma);

    delta_theta = abs(theta-theta0);   
    nu = delta_theta + nu0;
        
    %loop
    for i = 1:length(x)
      
        %calculate the Mach numbers after each pressure tap location from Prandtl-Meyer expansion
        M(i) = BisectionAlgorithm_PrandtlMeyer(nu(i), M_init(1), M_init(2), gamma, tol)
        
        %calculate pressure at each pressure tap location from Prandtl-Meyer from Mach number
        pi_p2(i) = CalcPressureExpansion(M2, M(i), gamma);
    end 
    
    if jj == 1
%         ptaps_p15_sea = p;
        pi_pinf_p15_sea = p2_pinf.*pi_p2;
    else
%         ptaps_m15_sea = p;
        pi_pinf_m15_sea = p2_pinf.*pi_p2;
    end
end

%cp
% cp_p15_sea = 2/(gamma*Minf^2)*(ptaps_p15_sea/pinf_p15 - 1);
% cp_m15_sea = 2/(gamma*Minf^2)*(ptaps_m15_sea/pinf_m15-1); 

cp_p15_sea_n = 2/(gamma*Minf^2)*(pi_pinf_p15_sea - 1);
cp_m15_sea_n = 2/(gamma*Minf^2)*(pi_pinf_m15_sea - 1); 


%% plot

figure
plot([1:0.01:20],CalculatePrandtlMayerAngle([1:0.01:20], gamma))
title('Prandtl-Meyer function')

figure
hold on
plot(x/c,cp_p15(1:end-1),'ok','DisplayName','lower surface: experiment')
plot(x/c,cp_p15_lin,'-.k','DisplayName','lower surface: linearized theory')
% plot(x/c,cp_p15_sea,'-.k','DisplayName','\alpha = +1.5 deg shock-expansion approx.')
plot(x/c,cp_p15_sea_n,'-k','DisplayName','lower surface: shock-expansion theory')
plot(x/c,cp_m15(1:end-1),'or','DisplayName','upper surface: experiment')
plot(x/c,cp_m15_lin,'-.r','DisplayName','upper surface: linearized theory')
% plot(x/c,cp_m15_sea,'-.r','DisplayName','\alpha = -1.5 deg shock-expansion approx.')
plot(x/c,cp_m15_sea_n,'-r','DisplayName','upper surface: shock-expansion theory')
title('pressure coefficient profiles')
xlim([0 1])
ylabel('cp [-]')
xlabel('x/c [-]')
legend('show','Location','northeast')
grid on



%% functions

function p = CalcStaticPressure(M, p0, gamma)
%{
    Args:
        M:      Mach number
        p0:     stagnation/total pressure
        gamma:  isentropic exponent
    Returns:
        p:     static pressure
    %}
    p = p0 / (1 + (gamma - 1)/2 * M^2)^(gamma/(gamma-1));
end

function p = CalcPressureNormalShock(M, pinf, gamma)
%{
    Args:
        M:      Mach number
        gamma:  isentropic exponent
    Returns:
        pratio: static pressure ratio across shock (normal shock rel.)
    %}
    p = pinf * (1 + (2*gamma)/(gamma+1) * (M^2-1));
end

function tan_theta = CalcTanTheta(M, beta, gamma)
    tan_theta = 2*cot(beta)*((M^2*sin(beta)^2-1)/(M^2*(gamma+cos(2*beta))+2));
end

% theta-beta-M relation
function beta = BisectionAlgorithm_ThetaBetaM(theta, M, beta_low, beta_up, gamma, tol)
    beta_mid = (beta_low + beta_up) / 2;                 % initial Mach number
    while abs(tan(theta) - CalcTanTheta(M, beta_mid, gamma)) > tol 
        beta_mid*180/pi;
        if CalcTanTheta(M, beta_mid, gamma) > tan(theta) % if RHS is smaller than LHS: Mach number is too high
            beta_up = beta_mid;
        else
            beta_low = beta_mid;                      % if RHS is greater than LHS: Mach number is too low
        end
        beta_mid = (beta_up + beta_low) / 2;
    end
    beta = beta_mid;
end

function nu = CalculatePrandtlMayerAngle(M, gamma)
    nu = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1)/(gamma + 1) .* (M.^2 - 1))) - atan(sqrt(M.^2 - 1));
end

% Prandtl-Meyer function
function M = BisectionAlgorithm_PrandtlMeyer(nu, M_low, M_up, gamma, tol)
    M_mid = (M_low + M_up) / 2;                 % initial Mach number
    while abs(nu - CalculatePrandtlMayerAngle(M_mid, gamma)) > tol 
        if CalculatePrandtlMayerAngle(M_mid, gamma) > nu % if RHS is smaller than LHS: Mach number is too high
            M_up = M_mid;
        else
            M_low = M_mid;                      % if RHS is greater than LHS: Mach number is too low
        end
        M_mid = (M_up + M_low) / 2
    end
    M = M_mid;
end

function p2_p1 = CalcPressureExpansion(M1, M2, gamma)
    % p1, M1: before expansion
    % p2, M2: after expansion
    p2_p1 = ((1 + (gamma-1)/2 * M1^2)/(1 + (gamma-1)/2 * M2^2))^(gamma/(gamma-1));
end

