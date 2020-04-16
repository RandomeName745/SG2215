% Define pipe properties

C_f = 0.005;                % friction coefficient
D = 0.0064;                 % inner pipe diameter
L = 2.67;                   % pipe length
A1 = 19.5E-3;               % pipe diameter at inlet of Venturi section
A2 = 10E-3;                 % pipe diameter at outlet of Venturi section

% Define fluid properties

gamma = 1.4;                % isentropic exponent of air
R = 287.058;                % gas constant of air

% Setup calculation of inlet Mach number with bisection algorithm

M_init = [0, 2];            % initial lower and upper limit of Mach number: M_init = [M_low, M_up];
tol = 1E-6;                 % tolerance for bisection algorithm

%% 4_1 Apply bisection algorithm
M_in = BisectionAlgorithm_41(L, C_f, D, gamma, M_init(1), M_init(2), tol);
disp(['Inlet Mach number for choked flow (M_2 = 1) at outlet of tube: ' num2str(M_in)]);

%% 4_2 Calculate the static pressure at tube inlet which is required to obtain M_2 = 1 (choked) at outlet of tube. This pressure equals stagnation pressure in stagnation chamber.
pstar = 100E3;  % static pressure at outlet (critical state)
p = CalcStaticPressureFromCritState(M_in, pstar, gamma);
disp(['Inlet pressure for choked flow (M_2 = 1) at outlet of tube: ' num2str(p/1E3) ' kPa']);

%% 4_3 Calculate local Mach numbers and static pressures with bisection algorithm
% p0 = 449.384E3;
x = [0.076, 0.76, 1.52, 2.28, 2.4, 2.53, 2.6, 2.63];
M = zeros(1, length(x));
p = zeros(1, length(x));
for i = 1:length(x)
    M(i) = BisectionAlgorithm_41(L - x(i), C_f, D, gamma, M_init(1), M_init(2), tol);
    p(i) = CalcStaticPressureFromCritState(M(i), pstar, gamma);
    p0(i) = CalcStagnationPressure(M(i), p(i), gamma);
end

%% Lab
str = fileread('lab_data1/Exp1.txt');
str = strrep( str, ' =', ',' );
C = textscan( str,'%s%f%f%f%f%f%f%f%f','headerlines',0,'delimiter',',');
data1 = cell2mat(C(2:size(C,2)));
data1 = data1*1E3;

str = fileread('lab_data1/Exp2.txt');
str = strrep( str, ' =', ',' );
C = textscan( str,'%s%f','headerlines',0,'delimiter',',');
data2 = cell2mat(C(2:size(C,2)));
data2 = data2*1E3;

str = fileread('lab_data1/Exp3.txt');
str = strrep( str, ' =', ',' );
C = textscan( str,'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',0,'delimiter',',');
data3 = cell2mat(C(2:size(C,2)));
data3 = data3*1E3;

%%%%%% 6_1a

T0 = 293.15;
for i = 1:size(data1, 2)
    % First assume n1 being at stagnation conditions (stagnation chamber)
    p01 = data1(1, i);
    p2 = data1(3, i);
    T1 = T0;
    for j = 1:20       
        M2 = CalcMach(p2, p01, gamma);
        T2 = CalcTemp(T0, M2, gamma);
        a2 = sqrt(gamma * R * T2);
        u2 = M2*a2;
        rho2 = data1(3, i)/(R*T2);
        m_dot(i, j) = rho2 * u2 * A2;

        M1 = M2 * (A2/A1) * (data1(3, i)/data1(2, i)) * sqrt(T1/T2);

        % Now correct for flow not being stagnant at n1

        T1 = CalcTemp(T0, M1, gamma);
        p01 = CalcStagnationPressure(M1, data1(2, i), gamma);
    end
        
end

figure(611)
plot(data1(12,:)/1E3,m_dot(:,end), 'Color', 'b', 'Linewidth', 1.5)
grid on
xl = xlabel('Back pressure $p_B$ in kPa');
yl = ylabel(['Mass flow rate ' '$\dot{m}$ in kg/s']);
set([xl, yl],'Interpreter','latex');

%%%%%% 6_1b

figure(612)    
hold on

for i = 1:size(data1,2)    
	scatter(x, data1(4:11,i)/data1(1,i), 'x', 'Linewidth', 1.5, 'DisplayName', ['$p_B$ = ' num2str(data1(12,i)/1E3) ' kPa']) 
end
lgnd = legend('-DynamicLegend');
set(lgnd, 'Location', 'southwest')
grid on
xlim([0, 2.67])
ylim([0, 1])
xl = xlabel('Downstream distance  $x$ from inltet in m');
yl = ylabel(['Normalized static pressure $\overline{p}$']);
set([xl, yl, lgnd],'Interpreter','latex');

%%%%%% 6_2a
figure(621)
hold on
scatter(x, data2(4:11)/data2(1),'x', 'b', 'Linewidth', 1.5, 'DisplayName', 'Experiment')
scatter(x, p/p0(1),'x', 'r', 'Linewidth', 1.5, 'DisplayName', 'Calculation')
grid on
lgnd = legend('-DynamicLegend');
xlim([0, 2.67])
ylim([0, 1])
xl = xlabel('Downstream distance  $x$ from inlet in m');
yl = ylabel(['Normalized static pressure $\overline{p}$']);
set([xl, yl, lgnd],'Interpreter','latex');


%%%%%% 6_2b
for i = 4:11
    M_62b(i-3) = BisectionAlgorithm_62(data2(i), pstar, gamma, M_init(1), M_init(2), tol);
end
figure(622)
hold on
scatter(x, M_62b, 'x', 'b', 'Linewidth', 1.5, 'DisplayName', 'Experiment')
scatter(x, M,'x', 'r', 'Linewidth', 1.5, 'DisplayName', 'Calculation')
grid on
lgnd = legend('-DynamicLegend');
xlim([0, 2.67])
ylim([0, 1])
xl = xlabel('Downstream distance  $x$ from inlet in m');
yl = ylabel(['Mach number $M$']);
set([xl, yl, lgnd],'Interpreter','latex');

%%%%%% 6_3c
for i = 1:size(data3,2)-1
    M_63c_isentropic(i) = CalcMach(pstar, data3(end,i), gamma);
    M_63c_rayleigh(i) = BisectionAlgorithm_63(data3(end,i), pstar, gamma, M_init(1), M_init(2), tol);    
end
figure(633)
hold on
plot(linspace(-D/2, D/2, size(data3,2)-1)*1E3, M_63c_isentropic, 'Color', 'b', 'Linewidth', 1.5, 'DisplayName', 'Isentropic relation')
plot(linspace(-D/2, D/2, size(data3,2)-1)*1E3, M_63c_rayleigh, 'Color', 'r', 'Linewidth', 1.5, 'DisplayName', 'Rayleigh Pitot tube formula')
lgnd = legend('-DynamicLegend');
grid on
xlim([-D/2, D/2]*1E3)
ylim([1.05, 1.25])
xl = xlabel('Radial position  $r$ at tube outlet in mm');
yl = ylabel(['Mach number $M$']);
set([xl, yl, lgnd],'Interpreter','latex');


%% Define functions
function M = BisectionAlgorithm_41(L, C_f, D, gamma, M_low, M_up, tol)
    M_mid = (M_up + M_low) / 2;                 % initial Mach number
    while abs(4*C_f*L/D - CalcRHS(M_mid, gamma)) > tol 
        M_mid = (M_up + M_low) / 2;
        if CalcRHS(M_mid, gamma) < 4*C_f*L/D    % if RHS is smaller than LHS: Mach number is too high
            M_up = M_mid;
        else
            M_low = M_mid;                      % if RHS is greater than LHS: Mach number is too low
        end
    end
    M = M_mid;
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