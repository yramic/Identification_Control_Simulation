%Author: Yannick Ramic
%Date: 27.07.2022


clear all;
clc;

%% Input Data:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Parameter');
addpath('Theta');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data_18 = importdata('Experiment_18.mat');
par.Data = Data_18;
% Aufbau Data_18:
% 1.Spalte: Time t
% 2.Spalte: T_HS_actual
% 3.Spalte: T_ICB_actual
% 4.Spalte: T_Water
% 5.Spalte: T_amb
% 6.Spalte: I_Peltier_actual
% 7.Spalte: s_Door (Door Status)
% 8.Spalte: s_Fan (Fan Status)
% 9.Spalte: I_ref_actual

t_ref = Data_18(:,1);

load('theta_a4_optimal_open.mat')
load('f_1_a4.mat')

par.theta_vektor = theta_vektor;
par.f_1 = f_1;

%% Parameter:

%Alle values for the System Coefficients need to be clarified in here!

par.m_air = 0.153; %kg
par.m_hs = 0.5; %kg
par.m_wall = 20; %kg
par.m_cargo = 1.7; %kg -> Aluminium: not insulated

par.c_air = 718; %J/(kg K)
par.c_hs = 910; %J/(kg K) -> Aluminium
par.c_wall = 1000; %J/(kg K)
par.c_cargo = 910; %J/(kg K) -> Aluminium

%Hinweis zu den Flächen: Bei der Wand wird sowohl aussen als auch innen die
%selbe Fläche angenommen, da der Unterschied marginal ist!
par.A_hs = 0.075; %m^2
par.A_wall = 1.5; %m^2
par.A_cargo = 0.0256; %kg -> A = 80x80 mm


par.lambda_wall = 0.15; %W/(m K)

par.d_wall = 0.03; %m

par.k_cargo = 120; %W/(m^2 K) -> Air / Aluminum
par.k_air_cooler = 120; %W/(m^2 K) -> Air / Aluminium (120)
par.k_air_wall = 15; %W/(m^2 K) -> Air / Wall
par.k_air_wall_i = 15; %W/(m^2 K) -> Wall inside

% par.k_air_cooler_erzwungen = 120;% 150.101; %W/(m^2 K) -> Air / Aluminium
% par.k_air_cooler_frei = 50; %W/(m^2 K) -> Air / Aluminium

par.k_air_cooler_erzwungen = 150.101; %W/(m^2 K) -> Air / Aluminium
par.k_air_cooler_frei = 50; %W/(m^2 K) -> Air / Aluminium

% par.T_atm = 20+273.15; %Degree Celsius
% par.T_water = 20+273.15; %Degree Celsius

% par.T_atm = 25+273.15; %Degree Celsius
% par.T_water = 15+273.15; %Degree Celsius

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis of T_Water:

Water_Vector = Data_18(:,4) + 273.15;
par.T_water_mean = mean(Water_Vector,'all');

% Analysis of T_amb:

Ambient_Temp_Vector = Data_18(:,5) + 273.15;
par.T_amb_mean = mean(Ambient_Temp_Vector,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.P_fan = 3; %W (0.5)
par.c_fan = 3; %W (0.5)

par.I = 4.7; %A (Amper) -> 4.7, 3.8
par.I_open = 3.8;

%Peltier Element:

par.U_max = 20; %V
par.I_max = 8.5; %A
par.delta_T_max = 71; %K or Degree Celsius

%First all necessary Parameters for the constats alpha, R & K need to get
%unpacked:

U_max = par.U_max;
I_max = par.I_max;
delta_T_max = par.delta_T_max;
T_water_mean = par.T_water_mean;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_hot = T_water + 273.15;
% par.T_hot = T_hot;
% 
% alpha_PE = U_max / T_hot;
% R_PE = ((T_hot - delta_T_max) * U_max) / (T_hot * I_max);
% K_PE = ((T_hot - delta_T_max) * U_max * I_max) / (2 * T_hot * delta_T_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha_PE = U_max / T_water_mean;
R_PE = ((T_water_mean - delta_T_max) * U_max) / (T_water_mean * I_max);
K_PE = ((T_water_mean - delta_T_max) * U_max * I_max) / (2 * T_water_mean * delta_T_max);

par.alpha_PE = alpha_PE;
par.R_PE = R_PE;
par.K_PE = K_PE;

%% Unpack all Parameters:

m_air = par.m_air; 
m_hs = par.m_hs; 
m_wall = par.m_wall; 
m_cargo = par.m_cargo;

c_air = par.c_air;
c_hs = par.c_hs; 
c_wall = par.c_wall; 
c_cargo = par.c_cargo;

A_hs = par.A_hs;
A_wall = par.A_wall;
A_cargo = par.A_cargo;


lambda_wall = par.lambda_wall;

d_wall = par.d_wall;

k_cargo = par.k_cargo;
k_air_cooler = par.k_air_cooler;
k_air_wall = par.k_air_wall;

T_amb_mean = par.T_amb_mean;

k_air_cooler_erzwungen = par.k_air_cooler_erzwungen;
k_air_cooler_frei = par.k_air_cooler_frei;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With Cargo : s_cargo = 1
% Without Cargo: s_cargo = 0
par.s_cargo = 0;
s_cargo = par.s_cargo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_fan = par.P_fan;
c_fan = par.c_fan;


I = par.I;





%% System Coefficients:

% Due to the fact that there are two cooling units, I adapted the
% parameters with the factor 2!


par.a1 = (A_wall*lambda_wall)/(m_wall*c_wall*d_wall);
par.a2 = (A_wall*k_air_wall)/(m_wall*c_wall);
par.a3 = (k_air_wall*A_wall)/(m_wall*c_wall);
par.a4 = (A_wall*lambda_wall)/(m_wall*c_wall*d_wall);
par.a5 = (A_cargo*k_cargo)/(m_cargo*c_cargo);

par.a6 = (2*k_air_cooler_erzwungen*A_hs)/(m_hs*c_hs);
par.a_frei_1 = (2*k_air_cooler_frei*A_hs)/(m_hs*c_hs);

% par.a6 = (2*k_air_cooler*A_hs)/(m_hs*c_hs);

par.a7 = (2*alpha_PE)/(m_hs*c_hs); %!!! Linearisierung
par.a8 = (2*K_PE)/(m_hs*c_hs);
% par.9 siehe unten
par.a10 = (A_cargo*k_cargo)/(m_air*c_air);
par.a9 = (A_wall*k_air_wall)/(m_air*c_air);
% par.a11 = (2*k_air_cooler*A_hs)/(m_air*c_air);
par.a11 = (2*k_air_cooler_erzwungen*A_hs)/(m_air*c_air);
par.a_frei_2 = (2*k_air_cooler_frei*A_hs)/(m_air*c_air);


% par.b1 = ...
par.b2 = (2*R_PE)/(2*m_hs*c_hs);
par.b3 = (2*c_fan)/(m_air*c_air);

par.e1 = -par.a3;
par.e2 = -par.a8;

%Unpack Parameter:
a1 = par.a1;
a2 = par.a2;
a3 = par.a3;
a4 = par.a4;
a5 = par.a5;
a6 = par.a6;
a7 = par.a7;
a8 = par.a8;
a9 = par.a9;
a10 = par.a10;
a11 = par.a11;

b2 = par.b2;
b3 = par.b3;

e1 = par.e1;
e2 = par.e2;

c(1) = a3;
c(2) = a4;
% c(2) = 8.346e-17;
par.c = c;


%% Find x1.stat:

% First it is necessary to analyze I whenever the Door is open:


% Assumptions to find a x1.stat:
% T_ICB_stat = 20 + 273.15;

T_ICB_stat = 16 + 273.15;
I_stat = 3.8;
%I_stat(T_ICB = 5) = 4.7 -> 1 Kühlaggregat
%I_stat(T_ICB = 5) = 3.8 -> 2 Kühlaggregat

%Formulate Equations:
syms T_wi_stat T_wo_stat T_cargo_stat T_HS_stat xi_stat
eqn_1 = a1*(T_wo_stat - T_wi_stat) - a2*(T_wi_stat - T_ICB_stat) == 0; 
eqn_2 = a3*(T_amb_mean - T_wo_stat) - a4*(T_wo_stat - T_wi_stat) == 0; 
eqn_3 = -a5*(T_cargo_stat - T_ICB_stat) == 0;
eqn_4 = a6*(T_ICB_stat - T_HS_stat) - a7*T_HS_stat*I_stat + ... 
    b2*(I_stat^2) + a8*(T_water_mean-T_HS_stat) == 0;
eqn_5 = xi_stat*abs(T_amb_mean - T_ICB_stat) + a9*(T_wi_stat - T_ICB_stat) + ...
    a10*(T_cargo_stat - T_ICB_stat) - a11*(T_ICB_stat - T_HS_stat) == 0;

% Simplify the system of Equations into a Matrix and a Vektor:
[Matrix_A, Vector_b] = equationsToMatrix([eqn_1, eqn_2, eqn_3, eqn_4, ...
    eqn_5], [T_wi_stat, T_wo_stat, T_cargo_stat, T_HS_stat, xi_stat]);

% Solve the system of Equations:
Solution_xi_stat_sym = linsolve(Matrix_A,Vector_b);

Solution_xi_stat = [];
Solution_xi_stat(5,1) = Solution_xi_stat_sym(5,1);
%Spalte 1 -> T_wi_stat
%Spalte 2 -> T_wo_stat
%Spalte 3 -> T_cargo_stat
%Spalte 4 -> T_HS_stat
%Spalte 5 -> xi_stat

for i = 1:4
    Solution_xi_stat(i,1) = Solution_xi_stat_sym(i,1)-273.15;  
end


par.g = 9.81; %Gravitation (m/s^2)
par.H = 0.5; %m
par.L = 0.5; %m
par.W = 0.5; %m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.W_AP = 0.5/2; %m (Assumption!!!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.c_D = 1; %Discharge Coefficient -> LaFayer Paper

par.c1_v = (par.g * par.H) / (8*par.L);
par.c2_v = - 9 / (8*par.L*(par.c_D ^2));
par.c3_v = (2*par.W_AP) / (par.L*par.W);
par.c4_v = 1 / (par.m_air*par.c_air);
par.c5_v = par.c3_v / par.c4_v;

% x1.stat:
x1_stat_cal = Solution_xi_stat(5,1) / (2*par.c3_v);
par.x1_stat = x1_stat_cal; %See LaFayer Paper (x1.stat = 0.1)
x1_stat = par.x1_stat;

par.xi = 2 * par.c3_v * par.x1_stat;
xi = par.xi;

% par.s_door = 0; %Tür geschlossen
% par.a9 = -par.xi;
% par.e3 = -par.a9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Initial Condition Problem + Numerical Coefficients:

% x0_x1 = 25 + 273.15;
% x0_x2 = 25 + 273.15;
% x0_x3 = 25 + 273.15;
% x0_x4 = 25 + 273.15;
% x0_x5 = 25 + 273.15;
% x0_x6 = 0;

% x1 ... T_w_i
% x2 ... T_w_o
% x3 ... T_cargo
% x4 ... T_HS
% x5 ... T_ICB
% x6 ... V/A

x0_x1 = Data_18(1,3) + 273.15; % Annahme: Start mit T_w_i = T_ICB
x0_x2 = 20 + 273.15; 
x0_x3 = 20 + 273.15;
x0_x4 = Data_18(1,2) + 273.15;
x0_x5 = Data_18(1,3) + 273.15;

x0 = [x0_x1; x0_x2; x0_x3; x0_x4; x0_x5];

%Time:

%All Time Variables are in seconds
t_beginn = Data_18(1,1);
t_end = Data_18(end,1);
tspan = [t_beginn,t_end];


%% ODE Solver:

options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',2.5); %'RelTol',1e-8
sol = ode15s(@(t,x) ICB(t,x,par), tspan, x0, options);

sim = (deval(sol,t_ref))';
sim(:,:) = sim(:,:)-273.15;
t_ref_h = t_ref./3600;

x_sim = (sol.y)';
t_sim = (sol.x)';
x_sim(:,1:5) = x_sim(:,1:5)-273.15;


%% Analyze I

% I_ref = Data_18(:,8);
% I = interp1(t_Data, I_ref, t, 'previous');
% 
% Interpolation = interp1(Data_18(:,1), Data_18(:,8), t_sim, 'previous');
% 
% figure
% plot(t_sim(:)./3600,Interpolation,'Linewidth',1)
% hold on
% plot(t_ref_h(:), Data_18(:,8))
% hold off

%% Linearisierung von Q_Punkt (Heatflow):

% NOTE! This section is activated when I want to have the stationary point
% for T_ICB_0 as well as for I_0:

T_HS_lin_1 = 0 + 273.15; % Linearization Point (Arbeitspunkt)
T_HS_lin_2 = T_HS_lin_1 + 10;
T_HS_lin_3 = T_HS_lin_1 - 10;

T_HS_lin_tot = [T_HS_lin_1, T_HS_lin_2, T_HS_lin_3];
I_lin_stat_tot = [];
Q_dot_lin_tot = [];
Q_dot_tot = [];
x_lin_tot = [];

I_lin = linspace(0,10,101);
    
th_1 = theta_vektor(1)*f_1(1);
th_2 = theta_vektor(2)*f_1(2);
th_3 = theta_vektor(3)*f_1(3);
th_4 = theta_vektor(4)*f_1(4);
th_5 = theta_vektor(5)*f_1(5);
    
th_6 = theta_vektor(6)*f_1(6);
th_7 = theta_vektor(7)*f_1(7);
th_8 = theta_vektor(8)*f_1(8);
    
th_9 = theta_vektor(9)*f_1(9);
th_10 = theta_vektor(10)*f_1(10);
th_11 = theta_vektor(11)*f_1(11);
th_12 = theta_vektor(12)*f_1(12);
th_13 = theta_vektor(13)*f_1(13);
    
par.th_1 = th_1;
par.th_2 = th_2;
par.th_3 = th_3;
par.th_4 = th_4;
par.th_5 = th_5;
par.th_6 = th_6;
par.th_7 = th_7;
par.th_8 = th_8;
par.th_9 = th_9;
par.th_10 = th_10;
par.th_11 = th_11;
par.th_12 = th_12;
par.th_13 = th_13;
    
%Formulate Equations:
s_fan_lin = 1; % Fan an im stationären Betrieb
s_door_lin = 0; % Tür zu
par.s_fan_lin = s_fan_lin;
par.s_door_lin = s_door_lin;

plot_stationary_HS = false;



if plot_stationary_HS == true

    for i = 1:numel(T_HS_lin_tot)
    
        par.T_HS_lin = T_HS_lin_tot(i);
        
        % fsolve is necessary to find I_stat for a HS Temperature of T_HS = -5
        % degree celsius:
        
        fun = @(x) eqn(x,par);
        x0_lin = [(0+273.15),(20+273.15),3,(0+273.15)];
        x_lin = fsolve(fun,x0_lin);
        I_lin_stat = x_lin(3); % If values are close use round !
    
        x_lin_tot(i,:) = x_lin;
        I_lin_stat_tot(i) = I_lin_stat;
        
        Q_dot_lin = (th_6*T_HS_lin_tot(i)-2*th_7*I_lin_stat)*I_lin + ...
            th_6*T_HS_lin_tot(i)*I_lin_stat - th_7*I_lin_stat^2 - ...
            (th_6*T_HS_lin_tot(i) - 2*th_7*I_lin_stat)*I_lin_stat;
        Q_dot = th_6*T_HS_lin_tot(i)*I_lin - ...
            th_7*(I_lin.^2);
    
        Q_dot_lin_tot(i,:) = Q_dot_lin;
        Q_dot_tot(i,:) = Q_dot;
    end
    
    % The Peltier Element is only used in between 3 - 7 A and at 0 but not in
    % between 0 - 3 A !
    find_3A = find(I_lin == 3);
    find_7A = find(I_lin == 7);
    
    % For a linearized Line trough zero I need 2 Points:
    linearization_point_1 = [0,0];
    half_Ampere = 3 + (7-3)/2;
    find_half = find(I_lin == half_Ampere);
    linearization_point_2 = [I_lin(find_half),Q_dot_lin_tot(1,find_half)];
    linearization_point_3 = [10, 2*Q_dot_lin_tot(1,find_half)];
    linearization_point_4 = [7, (7/half_Ampere)*Q_dot_lin_tot(1,find_half)];
    
    figure('Name','Heat Flow')
    p(1) = plot(I_lin,Q_dot_tot(1,:),'Color','k');
    hold on
    p(2) = plot(I_lin,Q_dot_tot(2,:),'Color','k');
    p(3) = plot(I_lin,Q_dot_tot(3,:),'Color','k');
    p(4) = plot(I_lin(find_3A:find_7A),Q_dot_lin_tot(1,find_3A:find_7A),'Color', ...
        'b','LineWidth',1.5,  'DisplayName', '$Q\dot{\,}$');
    p(5) = plot(I_lin(find_3A:find_7A),Q_dot_lin_tot(2,find_3A:find_7A),'Color', ...
        'b','LineWidth',1.5);
    p(6) = plot(I_lin(find_3A:find_7A),Q_dot_lin_tot(3,find_3A:find_7A),'Color', ...
        'b','LineWidth',1.5);
    p(7) = plot([linearization_point_1(1),linearization_point_4(1)], ...
        [linearization_point_1(2), linearization_point_4(2)], ...
        'LineWidth',2, 'Color', 'r',  'DisplayName', '$Q\dot{\,}_{lin.total}$');
    p(8) = xline(3,':');
    p(9) = xline(7,':');
    p(10) = plot(I_lin(1:find_3A),Q_dot_lin_tot(1,1:find_3A),'Color', ...
        'b','LineStyle', '--', 'DisplayName', '$Q\dot{\,}_{lin}$');
    p(11) = plot(I_lin(1:find_3A),Q_dot_lin_tot(2,1:find_3A),'Color', ...
        'b','LineStyle', '--');
    p(12) = plot(I_lin(1:find_3A),Q_dot_lin_tot(3,1:find_3A),'Color', ...
        'b','LineStyle','--');
    hold off
    plt_idx=[4,7,10];
    xlim([0,10])
    ylabel('Heat Flow $Q\dot{\,}$', 'Interpreter', 'latex')
    xlabel('Electrical Current I', 'Interpreter', 'latex')
    legend([p(plt_idx)],{'$Q\dot{\,}$', '$Q\dot{\,}_{lin.total}$', '$Q\dot{\,}_{lin}$'}, ...
        'Interpreter', 'latex')
    
    % Now I want to have the Parameters for the linearized Q_dot:
    lambda_Q_dot = 2*Q_dot_lin_tot(1,find_half) / 10;
    par.lambda_Q_dot = lambda_Q_dot;
end

%% Save Parameters:

% First T_wall_inside:
T_wall_inside = sim(:,1); %In Degree Celsius

T_wall_inside_1(:,1) = T_wall_inside(1:3131,1);
T_wall_inside_2(:,1) = T_wall_inside(3132:(3131*2),1);
T_wall_inside_3(:,1) = T_wall_inside(((3131*2)+1):end,1);

% save('T_wall_inside.mat','T_wall_inside');
% save('T_wall_inside_1.mat','T_wall_inside_1');
% save('T_wall_inside_2.mat','T_wall_inside_2');
% save('T_wall_inside_3.mat','T_wall_inside_3');

% T_wall_outside:

T_wall_outside = sim(:,2); %In Degree Celsius

T_wall_outside_1(:,1) = T_wall_outside(1:3131,1);
T_wall_outside_2(:,1) = T_wall_outside(3132:(3131*2),1);
T_wall_outside_3(:,1) = T_wall_outside(((3131*2)+1):end,1);

%This is necessary for the Identification!
% save('T_wall_outside.mat','T_wall_outside');
% save('T_wall_outside_1.mat','T_wall_outside_1');
% save('T_wall_outside_2.mat','T_wall_outside_2');
% save('T_wall_outside_3.mat','T_wall_outside_3');

% save('PARAMETER.mat','par');


%% Additional analysis!
% Calculate the R-squared value for each column of data
data_measured = Data_18(:,3);
data_result = sim(:,5);
R2 = zeros(1, size(data_measured, 2));
for i = 1:size(data_measured, 2)
    y = data_measured(:, i);
    ypred = data_result(:, i);
    SSres = sum((y - ypred).^2);
    SStot = sum((y - mean(y)).^2);
    R2(i) = 1 - SSres/SStot;
end

% Calculate the overall R-squared value as the average of the individual R-squared values
overall_R2 = mean(R2);

% Display the individual and overall R-squared values
disp(['Individual R-squared values: ', num2str(R2)])
disp(['Overall R-squared value: ', num2str(overall_R2)])

% Calculate NRMSE
range = max(data_measured) - min(data_measured);
nrmse = sqrt(mean((data_result - data_measured).^2)) / range;
fit_nrmse = 1-nrmse;
disp(['Overall NRMSE value: ', num2str(overall_R2)])

%% Plots

% Zusammenfassung der Simulation MIT DEVAL!:
figure('Name','Simulation with Kelvin')
plot(Data_18(:,1)./3600,Data_18(:,5),'Linewidth',1)
hold on
plot(Data_18(:,1)./3600,Data_18(:,4),'Linewidth',1)
plot(t_ref_h(:),sim(:,1),'Linewidth',1)
plot(t_ref_h(:),sim(:,2),'Linewidth',1)
plot(t_ref_h(:),sim(:,3),'Linewidth',1)
plot(t_ref_h(:),sim(:,4),'Linewidth',1)
plot(t_ref_h(:),sim(:,5),'Linewidth',1)
hold off
xlim([t_beginn,t_end/3600])
ylim([round(min(sim(:,1:5),[],'all')-1,0),max(sim(:,1:5),[],'all')+5])
% xticks([0, 0.5*3600, 3600, 1.5*3600, 2*3600, 2.5*3600, 3*3600, 3.5*3600,...
%     4*3600, 4.5*3600, 5*3600, 5.5*3600, 6*3600, 6.5*3600, 7*3600, ...
%     7.5*3600, 8*3600, 8.5*3600, 9*3600, 9.5*3600, 10*3600])
% xticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5',...
%     '5.5','6','6.5', '7', '7.5', '8', '8.5', '9', '9.5', '10'})
xlabel('Time [h]')
ylabel('Temperature [\circC]')
legend('T_{amb}','T_{Water}' ,'T_{wall.inside}','T_{wall.outside}','T_{Cargo}','T_{HS}',...
    'T_{ICB}')
title('Simulation with Kelvin')


% Vergleich der ICB Temperatur:
figure('Name','Difference T_ICB')
plot(t_ref_h(:),sim(:,5),'Linewidth',2)
hold on
plot(t_ref_h(:),Data_18(:,3),'Linewidth',1)
hold off
xlim([t_beginn,t_end/3600])
ylim([round(min(Data_18(:,2),[],'all')-1,0),max(sim(:,4),[],'all')+5])
xlabel('Time [h]')
ylabel('Temperature [\circC]')
legend('T_{Simulation}','T_{Actual}')
title('Difference T_{ICB}')

% Vergleich der HS Temperatur:
figure('Name','Difference T_HS')
plot(t_ref_h(:),sim(:,4),'Linewidth',2)
hold on
plot(t_ref_h(:),Data_18(:,2),'Linewidth',1)
hold off
xlim([t_beginn,t_end/3600])
ylim([round(min(sim(:,4),[],'all')-1,0),max(sim(:,4),[],'all')+5])
xlabel('Time [h]')
ylabel('Temperature [\circC]')
legend('T_{Simulation}','T_{Actual}')
title('Difference T_{HS}')

%% Functions:

function dxdt = ICB(t,x,par)


%unpack all Parameters:
a1 = par.a1;
a2 = par.a2;
a3 = par.a3;
a4 = par.a4;
a5 = par.a5;
a6 = par.a6;
a7 = par.a7;
a8 = par.a8;
a9 = par.a9;
a10 = par.a10;
a11 = par.a11;

a_frei_1 = par.a_frei_1;
a_frei_2 = par.a_frei_2;


% b1 = par.b1;
b2 = par.b2;
b3 = par.b3;


% T_atm = par.T_atm;
% T_water = par.T_water;

s_cargo = par.s_cargo;
c = par.c;
f_1 = par.f_1;
th = par.theta_vektor;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data_18 = par.Data;
t_Data = Data_18(:,1);

s_Door_ref = Data_18(:,7);
s_Door = interp1(t_Data, s_Door_ref, t, 'previous');

I_ref = Data_18(:,6);
I = interp1(t_Data, I_ref, t, 'previous');

s_Fan_ref = Data_18(:,8);
s_Fan = interp1(t_Data, s_Fan_ref, t, 'previous');

T_water_ref = Data_18(:,4)+273.15;
T_water = interp1(t_Data, T_water_ref, t, 'previous');

T_amb_ref = Data_18(:,5)+273.15;
T_amb = interp1(t_Data, T_amb_ref, t, 'previous');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1_v = par.c1_v;
c2_v = par.c2_v;
c3_v = par.c3_v;
c4_v = par.c4_v;
x1_stat = par.x1_stat;
xi = par.xi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T_water = par.T_water;

%Beschreibung des Vektors x:
% x(1) -> T_wall_inside
% x(2) -> T_wall_outside
% x(3) -> T_cargo
% x(4) -> T_hs
% x(5) -> T_ICB
% x(6) -> (V_punkt / A_ap) [m/s]


%All 6 Differential Equations
dxdt(1) = th(1)*f_1(1)*(x(2)-x(1)) - th(2)*f_1(2)*(x(1)-x(5));
dxdt(2) = th(3)*f_1(3)*(T_amb-x(2)) - c(2)*(x(2)-x(1));
dxdt(3) = 0;
dxdt(4) = (th(4)*f_1(4)*s_Fan + (1-s_Fan)*th(5)*f_1(5))*(x(5)-x(4)) -...
    th(6)*f_1(6)*I*x(4) + ...
    th(7)*f_1(7)*(I^2) +...
    th(8)*f_1(8)*(T_water-x(4));
dxdt(5) = th(9)*f_1(9)*(x(1)-x(5)) - ...
    (th(10)*f_1(10)*s_Fan + (1-s_Fan)*th(11)*f_1(11))*(x(5)-x(4)) + ...
    th(12)*f_1(12)*(T_amb-x(5))*s_Door +...
    th(13)*f_1(13)*s_Fan;

dxdt = dxdt';
end

% Necessary Function to find a stationary value for the current I for a
% given HS Temperature:

function F = eqn(x, par)
    s_fan_lin = par.s_fan_lin;
    s_door_lin = par.s_door_lin;
    T_HS_lin = par.T_HS_lin;
    T_water_mean = par.T_water_mean;
    T_amb_mean = par.T_amb_mean;
    c = par.c;

    th_1 = par.th_1;
    th_2 = par.th_2;
    th_3 = par.th_3;
    th_4 = par.th_4;
    th_5 = par.th_5;
    th_6 = par.th_6;
    th_7 = par.th_7;
    th_8 = par.th_8;
    th_9 = par.th_9;
    th_10 = par.th_10;
    th_11 = par.th_11;
    th_12 = par.th_12;
    th_13 = par.th_13;

    % x(1) ... T_wi
    % x(2) ... T_wo
    % x(3) ... I
    % x(4) ... T_ICB

    F(1) = th_1*(x(2) - x(1)) - th_2*(x(1) - x(4));
    F(2) = th_3*(T_amb_mean - x(2)) - c(2)*(x(2) - x(1)); 
    F(3) = (th_4*s_fan_lin + (1-s_fan_lin)*th_5)*(x(4) - T_HS_lin) - ...
        th_6*T_HS_lin*x(3) + ... 
        th_7*(x(3)^2) + th_8*(T_water_mean-T_HS_lin);
    F(4) = th_9*(x(1) - x(4)) - ...
        (th_10*s_fan_lin + (1-s_fan_lin)*th_11)*(x(4) - T_HS_lin) + ...
        th_12*(T_amb_mean - x(4))*s_door_lin + ...
        th_13*s_fan_lin;
end

% Necessary Function to find a stationary value for the current I for a
% given ICB Temperature, since this is the value that needs to be 
% controlled:

function F_ICB = eqn_ICB(x, par)
    s_fan_lin = par.s_fan_lin;
    s_door_lin = par.s_door_lin;
    T_ICB_lin = par.T_ICB_lin;
    T_water_mean = par.T_water_mean;
    T_amb_mean = par.T_amb_mean;
    c = par.c;

    th_1 = par.th_1;
    th_2 = par.th_2;
    th_3 = par.th_3;
    th_4 = par.th_4;
    th_5 = par.th_5;
    th_6 = par.th_6;
    th_7 = par.th_7;
    th_8 = par.th_8;
    th_9 = par.th_9;
    th_10 = par.th_10;
    th_11 = par.th_11;
    th_12 = par.th_12;
    th_13 = par.th_13;

    % x(1) ... T_wi
    % x(2) ... T_wo
    % x(3) ... I
    % x(4) ... T_HS

    F_ICB(1) = th_1*(x(2) - x(1)) - th_2*(x(1) - T_ICB_lin);
    F_ICB(2) = th_3*(T_amb_mean - x(2)) - c(2)*(x(2) - x(1)); 
    F_ICB(3) = (th_4*s_fan_lin + (1-s_fan_lin)*th_5)*(T_ICB_lin - x(4)) - ...
        th_6*x(3)*x(4) + ... 
        th_7*(x(3)^2) + th_8*(T_water_mean-x(4));
    F_ICB(4) = th_9*(x(1) - T_ICB_lin) - ...
        (th_10*s_fan_lin + (1-s_fan_lin)*th_11)*(T_ICB_lin - x(4)) + ...
        th_12*(T_amb_mean - T_ICB_lin)*s_door_lin + ...
        th_13*s_fan_lin;
end
