% Constants
Minf = 2;
pinf = 1*10^5;
AOA = 1.5;              % angle of attack
gamma = 1.4;

%%

% Wing properties
r = 210;                % radius for wing curvature
lu_deg = 360/(2*pi*r);  % degree (°) per length unit (l.u.)

x = linspace(37,-37);   % x-positions of the wing in l.u. (x(1) = wing tip, x(end) = trailing edge)
phi = x.*lu_deg;        % x-positions of the wing transformed into

theta_u = phi - AOA;
theta_l = phi + AOA;

%%

% pinf = CalcStaticPressure(Minf,  p0, gamma);

cp_u = (2*theta_u) ./ (sqrt(Minf^2 - 1));
cp_l = (2*theta_l) ./ (sqrt(Minf^2 - 1));

LocalPressure(cp_u, gamma, Minf, pinf)

M_2_u = BisectionAlgorithm_M_2(theta_u, beta_u, M_init(1), M_init(2), gamma1, tol);
M_2_l = BisectionAlgorithm_M_2(theta_l, beta_l, M_init(1), M_init(2), gamma1, tol);
disp(['Experiment3: Upper Mach number behind shock (for beta_u = ' num2str(beta_u*180/pi) '° and theta_u = ' num2str(theta_u*180/pi) '°) : ' num2str(M_2_u)]);
disp(['Experiment3: Lower Mach number behind shock (for beta_l = ' num2str(beta_l*180/pi) '° and theta_u = ' num2str(theta_l*180/pi) '°) : ' num2str(M_2_l)]);


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

function tan_theta = CalcTanTheta(M, beta, gamma)
    tan_theta = 2*cot(beta)*((M^2*sin(beta)^2-1)/(M^2*(gamma+cos(2*beta))+2));
end

% theta-beta-M relation
function M_2 = BisectionAlgorithm_ThetaBetaM(theta, M, beta_low, beta_up, gamma, tol)
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
