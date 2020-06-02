% Constants
Minf = 2;
pinf = 1*10^5;
AOA = 1.5;              % angle of attack
gamma = 1.4;

%%

% Wing properties
r = 210;                % radius for wing curvature
lu_deg = 360/(2*pi*r);  % degree (�) per length unit (l.u.)

x = linspace(37,-37);   % x-positions of the wing in l.u. (x(1) = wing tip, x(end) = trailing edge)
phi = x.*lu_deg;        % x-positions of the wing transformed into

theta_u = phi - AOA;
theta_l = phi + AOA;

%%

% pinf = CalcStaticPressure(Minf,  p0, gamma);

cp_u = (2*theta_u) ./ (sqrt(Minf^2 - 1));
cp_l = (2*theta_l) ./ (sqrt(Minf^2 - 1));

LocalPressure(cp_u, gamma, Minf, pinf)

beta_init = [0, pi/2];
tol = 1E-6;

% Calculate shock angle for oblique shock at leading edge
beta_u = BisectionAlgorithm_ThetaBetaM(theta_u(1), Minf, beta_init(1), beta_init(2), gamma, tol);
beta_l = BisectionAlgorithm_ThetaBetaM(theta_l(1), Minf, beta_init(1), beta_init(2), gamma, tol);
% Calculate normal Mach number
Mn1_u = Minf*beta_u;
Mn1_l = Minf*beta_l;
% Calculate pressure ratio across oblqiue shock from normal shock relation
pratio_u =  CalcPressureRatioNormalShock(Mn1_u,gamma);
pratio_l =  CalcPressureRatioNormalShock(Mn1_l,gamma);
% Calculate normal Mach number behind oblique shock
Mn2_u = (1 + (gamma - 1)/2 * Mn1_u^2) / (gamma*Mn1_u^2 - (gamma - 1)/2);
Mn2_l = (1 + (gamma - 1)/2 * Mn1_l^2) / (gamma*Mn1_l^2 - (gamma - 1)/2);
% Calculate Mach number behind oblique shock
M2_u = Mn2_u / sin(beta_u - theta_u(1));
M2_l = Mn2_l / sin(beta_l - theta_l(1));
% 
nu2_u = zeros(size(x));
nu2_l = zeros(size(x));
nu2_u(1) = CalculatePrandtlMayerAngle(M2_u, gamma);
nu2_l(1) = CalculatePrandtlMayerAngle(M2_l, gamma);




function p = LocalPressure(cp, gamma, Minf, pinf)
    p = 0.5*(gamma*Minf^2*pinf*cp) + pinf;
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

function pratio = CalcPressureRatioNormalShock(M, gamma)
%{
    Args:
        M:      Mach number
        gamma:  isentropic exponent
    Returns:
        pratio: static pressure ratio across shock (normal shock rel.)
    %}
    pratio = 1 + (2*gamma)/(gamma+1) * (M^2-1);
end

function tan_theta = CalcTanTheta(M, beta, gamma)
    tan_theta = 2*cot(beta)*((M^2*sin(beta)^2-1)/(M^2*(gamma+cos(2*beta))+2));
end

% theta-beta-M relation
function beta = BisectionAlgorithm_ThetaBetaM(theta, M, beta_low, beta_up, gamma, tol)
    theta = theta*pi/180;
    beta_mid = (beta_low + beta_up) / 2;                 % initial Mach number
    while abs(tan(theta) - CalcTanTheta(M, beta_mid, gamma)) > tol 
        beta_mid = (beta_up + beta_low) / 2;
        if CalcTanTheta(M, beta_mid, gamma) > tan(theta) % if RHS is smaller than LHS: Mach number is too high
            beta_up = beta_mid;
        else
            beta_low = beta_mid;                      % if RHS is greater than LHS: Mach number is too low
        end
    end
    beta = beta_mid;
end

function nu = CalculatePrandtlMayerAngle(M, gamma)
    nu = sqrt((gamma + 1) / (gamma - 1)) * atan((gamma + 1)/(gamma - 1) * (M^2 - 1)) - atan(M^2 - 1);
end
