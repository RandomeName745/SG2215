% Define shock tube properties

% Define fluid properties

gamma1 = 1.4;               % isentropic exponent of air
R1 = 287.058;               % gas constant of air
T1 = 300;

gamma4 = 1.667;             % isentropic exponent of helium
R4 = 2076.9;                % gas constant of helium
T4 = 300;

tol = 1E-6;                 % tolerance for bisection algorithm


%% 5.1
M2 = 0.46;
p1 = 49*10^3;   %[Pa]
p2 = 96*10^3;   %[Pa]

% Setup calculation of Mach number Ms of shock with bisection algorithm
% Depends on required Mach number M2 behind shock
M_init = [0, 2];            % initial lower and upper limit of Mach number: M_init = [M_low, M_up];

% 1a
M_s = BisectionAlgorithm_M_s(M2, M_init(1), M_init(2), gamma1, tol);
disp(['Experiment1: Mach number of shock wave (for M2 = ' num2str(M2) ') : ' num2str(M_s)]);

% 1b
a1 = CalcSpeedOfSound(gamma1, R1, T1);
a4 = CalcSpeedOfSound(gamma4, R4, T4);
p4_p1_th = Calcp4_p1(gamma1, gamma4, M_s, a1, a4);
disp(['Experiment1: theoretically required pressure ratio p4/p1: ' num2str(p4_p1_th)]);

% 1c
p4_p1_calib = Calcp4_p1_calib(M_s);
disp(['Experiment1: actual required pressure ratio p4/p1: ' num2str(p4_p1_calib)]);

% 2a
M_s_exp = CalcMs_p2_p1(gamma1,p1,p2);
disp(['Experiment1: Mach number of shock wave (for p2/p1 = ' num2str(p2/p1) ') : ' num2str(M_s_exp)]);

% 2b
M_2_exp = CalcM2(M_s_exp, gamma1);
disp(['Experiment1: Mach number behind shock (for M_s_exp = ' num2str(M_s_exp) ') : ' num2str(M_2_exp)]);
disp(['Experiment1: Relative difference to target M2 : ' num2str(abs((1-M_2_exp/M2)*100)) ' % ']);

% 2c
% Observation referring to the obtained M2 value:
% No bow shock or oblique shock is formed at the wedge as M2<1
% Reflection of incident shockwave from wall can be observed when
% travelling back.

%% 5.2
M2 = 1.04;
p1 = 10*10^3;   %[Pa]
p2 = 52*10^3;   %[Pa]

% Setup calculation of Mach number Ms of shock with bisection algorithm
% Depends on required Mach number M2 behind shock
M_init = [0, 3];            % initial lower and upper limit of Mach number: M_init = [M_low, M_up];

% 1a
M_s = BisectionAlgorithm_M_s(M2, M_init(1), M_init(2), gamma1, tol);
disp(['Experiment2: Mach number of shock wave (for M2 = ' num2str(M2) ') : ' num2str(M_s)]);

% 1b
a1 = CalcSpeedOfSound(gamma1, R1, T1);
a4 = CalcSpeedOfSound(gamma4, R4, T4);
p4_p1_th = Calcp4_p1(gamma1, gamma4, M_s, a1, a4);
disp(['Experiment2: theoretically required pressure ratio p4/p1: ' num2str(p4_p1_th)]);

% 1c
p4_p1_calib = Calcp4_p1_calib(M_s);
disp(['Experiment2: actual required pressure ratio p4/p1: ' num2str(p4_p1_calib)]);

% 2a
M_s_exp = CalcMs_p2_p1(gamma1,p1,p2);
disp(['Experiment2: Mach number of shock wave (for p2/p1 = ' num2str(p2/p1) ') : ' num2str(M_s_exp)]);

% 2b
M_2_exp = CalcM2(M_s_exp, gamma1);
disp(['Experiment2: Mach number behind shock (for M_s_exp = ' num2str(M_s_exp) ') : ' num2str(M_2_exp)]);
disp(['Experiment2: Relative difference to target M2 : ' num2str(abs((1-M_2_exp/M2)*100)) ' % ']);

% 2c
% Observation referring to calculated M2 and theta-beta-M-curve
% M2 > 1 but theta > theta_max, thus bow shock which is detached from the wedge is formed behind
% the incident shock wave at the wedge
% Influence of theta angle can be seen as bow is formed accordingly (larger
% angle in upper region)
% theta > theta_max, thus solution can not be seen in theta-beta-M curve 
% (oblique shock only for theta <= theta_max)

%% 5.3
M2 = 1.34;
p1 = 10*10^3;   %[Pa]
p2 = 100*10^3;   %[Pa]

% from snapshot
d_thetabeta_u = 59.3;  %[deg]
d_thetabeta_l = 50.9;  %[deg]

% Setup calculation of Mach number Ms of shock with bisection algorithm
% Depends on required Mach number M2 behind shock
M_init = [0, 4];            % initial lower and upper limit of Mach number: M_init = [M_low, M_up];

% 1a
M_s = BisectionAlgorithm_M_s(M2, M_init(1), M_init(2), gamma1, tol);
disp(['Experiment3: Mach number of shock wave (for M2 = ' num2str(M2) ') : ' num2str(M_s)]);

% 1b
a1 = CalcSpeedOfSound(gamma1, R1, T1);
a4 = CalcSpeedOfSound(gamma4, R4, T4);
p4_p1_th = Calcp4_p1(gamma1, gamma4, M_s, a1, a4);
disp(['Experiment3: theoretically required pressure ratio p4/p1: ' num2str(p4_p1_th)]);

% 1c
p4_p1_calib = Calcp4_p1_calib(M_s);
disp(['Experiment3: actual required pressure ratio p4/p1: ' num2str(p4_p1_calib)]);


% 2a
M_s_exp = CalcMs_p2_p1(gamma1,p1,p2);
disp(['Experiment3: Mach number of shock wave (for p2/p1 = ' num2str(p2/p1) ') : ' num2str(M_s_exp)]);

% 2b
M_2_exp = CalcM2(M_s_exp, gamma1);
disp(['Experiment3: Mach number behind shock (for M_s_exp = ' num2str(M_s_exp) ') : ' num2str(M_2_exp)]);
disp(['Experiment3: Relative difference to target M2 : ' num2str(abs((1-M_2_exp/M2)*100)) ' % ']);

% 2c
theta_u = 8*pi/180; %[rad]
beta_u = d_thetabeta_u*pi/180  + theta_u;  %[rad]
theta_l = 4*pi/180; %[rad]
beta_l = d_thetabeta_l*pi/180 + theta_l;  %[rad]

M_2_u = BisectionAlgorithm_M_2(theta_u, beta_u, M_init(1), M_init(2), gamma1, tol);
M_2_l = BisectionAlgorithm_M_2(theta_l, beta_l, M_init(1), M_init(2), gamma1, tol);
disp(['Experiment3: Upper Mach number behind shock (for beta_u = ' num2str(beta_u*180/pi) '° and theta_u = ' num2str(theta_u*180/pi) '°) : ' num2str(M_2_u)]);
disp(['Experiment3: Lower Mach number behind shock (for beta_l = ' num2str(beta_l*180/pi) '° and theta_u = ' num2str(theta_l*180/pi) '°) : ' num2str(M_2_l)]);

% 2d
% Observation referring to calculated M2 and theta-beta-M-curve
% M2 > 1 but theta < theta_max, thus oblique shocks formated at the wedge tip 
% when the incident shock wave hits the wedge
% Influence of theta angle can be seen as bow is formed accordingly (larger beta for larger theta
% in upper part)
% theta < theta_max, thus solution can be seen in theta-beta-M curve: below
% M2=1 curve -> weak shock

% 2e
% Which is the most reliable value?
% Not the one from 2c because of inaccuracy when measuring beta


%% Define functions

% Shock Mach number
function M_s = BisectionAlgorithm_M_s(M2, M_low, M_up, gamma, tol)
    M_mid = (M_low + M_up) / 2;                 % initial Mach number
    while abs(M2 - CalcM2(M_mid, gamma)) > tol 
        M_mid = (M_up + M_low) / 2;
        if CalcM2(M_mid, gamma) > M2    % if RHS is smaller than LHS: Mach number is too high
            M_up = M_mid;
        else
            M_low = M_mid;                      % if RHS is greater than LHS: Mach number is too low
        end
    end
    M_s = M_mid;
end

function M_s = CalcMs_p2_p1(gamma1,p1,p2)

    M_s = sqrt((gamma1+1)/(2*gamma1)*(p2/p1-1)+1);
end

% Speed of sound
function a = CalcSpeedOfSound(gamma, R, T)
    a = sqrt(gamma*R*T);
end

% theta-beta-M relation
function M_2 = BisectionAlgorithm_M_2(theta, beta, M_low, M_up, gamma, tol)
    M_mid = (M_low + M_up) / 2;                 % initial Mach number
    while abs(tan(theta) - CalcTanTheta(M_mid, beta, gamma)) > tol 
        M_mid = (M_up + M_low) / 2;
        if CalcTanTheta(M_mid, beta, gamma) > tan(theta) % if RHS is smaller than LHS: Mach number is too high
            M_up = M_mid;
        else
            M_low = M_mid;                      % if RHS is greater than LHS: Mach number is too low
        end
    end
    M_2 = M_mid;
end


function p4_p1 = Calcp4_p1(gamma1, gamma4, M_s, a1, a4)
    p4_p1 = (2*gamma1*M_s^2 - (gamma1 - 1)) / (gamma1 + 1) * (1 - (gamma4 - 1) / (gamma1 + 1) * a1/a4 * (M_s - 1/M_s))^(-2*gamma4 / (gamma4 - 1));
end
function p4_p1_calib = Calcp4_p1_calib(M_s)
    p4_p1_calib = exp(0.31*M_s^3 - 2.6*M_s^2 + 8.1*M_s - 5.5);
end

function M2 = CalcM2(M_s, gamma)
    M2 = (2*(M_s^2 - 1)) / ((2*gamma*M_s^2 - (gamma - 1))^(1/2) *((gamma - 1)*M_s^2 + 2)^(1/2));
end

function tan_theta = CalcTanTheta(M, beta, gamma)
    tan_theta = 2*cot(beta)*((M^2*sin(beta)^2-1)/(M^2*(gamma+cos(2*beta))+2));
end

function M = BisectionAlgorithm_62(p, pstar, gamma, M_low, M_up, tol)
    M_mid = (M_up + M_low) / 2;                 % initial Mach number
    while abs(p/pstar - CalcStaticPressureFromCritState(M_mid, pstar, gamma)/pstar) > tol 
        M_mid = (M_up + M_low) / 2;
        if CalcStaticPressureFromCritState(M_mid, pstar, gamma)/pstar < p/pstar    % if RHS is smaller than LHS: Mach number is too high
            M_up = M_mid;
        else
            M_low = M_mid;                      % if RHS is greater than LHS: Mach number is too low
        end
    end
    M = M_mid;
end

function M = BisectionAlgorithm_63(p0, pstar, gamma, M_low, M_up, tol)
    M_mid = (M_up + M_low) / 2;                 % initial Mach number
    while abs(p0/pstar - CalcRHSRayleigh(M_mid, gamma)) > tol 
        M_mid = (M_up + M_low) / 2;
        if CalcRHSRayleigh(M_mid, gamma) > p0/pstar    % if RHS is smaller than LHS: Mach number is too high
            M_up = M_mid;
        else
            M_low = M_mid;                      % if RHS is greater than LHS: Mach number is too low
        end
    end
    M = M_mid;
end

function RHS = CalcRHS(M_in, gamma)
    %{
    Args:
        M_in:   inlet Mach number
        gamma:  isentropic exponent
    %}
    RHS = (1 - M_in^2)/(gamma*M_in^2) + (gamma + 1)/(2*gamma)*log(((gamma + 1)*M_in^2)/(2 + (gamma - 1)*M_in^2));
end

function RHS = CalcRHSRayleigh(M, gamma)
    %{
    Args:
        M:      Mach number
        gamma:  isentropic exponent
    %}
    RHS = (((gamma + 1)^2 * M^2) / (4*gamma*M^2-2*(gamma - 1)))^(gamma / (gamma - 1)) * ((1 - gamma + 2*gamma*M^2) / (gamma + 1));
end

function M = CalcMach(p, p0, gamma)
%{
    Args:
        p:     static pressure
        p0:     stagnation/total pressure
        gamma:  isentropic exponent
    Returns:
        M:      Mach number
    %}
    M = sqrt(((p0/p)^((gamma -1 )/gamma) -1 ) * 2 / (gamma -1));
end

function T = CalcTemp(T0, M, gamma)
%{
    Args:
        T0:     stagnation temperature
        M:      Mach number
        gamma:  isentropic exponent
    Returns:
        T:      Static temperature
    %}
    T = T0 / (1 + (gamma -1)/2 * M^2);
end

function p = CalcStaticPressureFromCritState(M, pstar, gamma)
    %{
    Args:
        M:      Mach number
        pstar:  critical state static pressure
        gamma:  isentropic exponent
    Returns:
        p:     	static pressure
    %}
    p = pstar/M * sqrt((gamma + 1)/(2 + (gamma -1)*M^2));
end

function p0 = CalcStagnationPressure(M, p, gamma)
%{
    Args:
        M:      Mach number
        p:      static pressure
        gamma:  isentropic exponent
    Returns:
        p0:     stagnation/total pressure
    %}
    p0 = p * (1 + (gamma - 1)/2 * M^2)^(gamma/(gamma-1));
end

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

function M1 = CalcMach1NormalShock(M2, gamma)
    M1 = sqrt( (1 + 0.5 * M2^2 * (gamma - 1)) / (gamma *M2^2 - 0.5 * (gamma - 1)) );
end