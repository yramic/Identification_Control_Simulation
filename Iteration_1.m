%Author: Yannick Ramic
%Date: 27.07.2022


clear all;
clc;

%% Input Data:

addpath('Theta');
addpath('Parameter_BSC');
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
%% Parameter:

%Alle values for the System Coefficients need to be clarified in here!

par.m_air = 0.153; %kg
par.m_hs = 0.5; %kg %Original: 0.5
par.m_wall = 20; %kg %Originial: 20
par.m_cargo = 1.7; %kg -> Aluminium: not insulated

par.c_air = 718; %J/(kg K)
par.c_hs = 1000; %J/(kg K) -> Aluminium
par.c_wall = 1000; %J/(kg K) %Original: 1000
par.c_cargo = 910; %J/(kg K) -> Aluminium

%Hinweis zu den Flächen: Bei der Wand wird sowohl aussen als auch innen die
%selbe Fläche angenommen, da der Unterschied marginal ist!
par.A_hs = 0.075; %m^2
par.A_wall = 1.5; %m^2
par.A_cargo = 0.0256; %kg -> A = 80x80 mm


par.lambda_wall = 0.15; %W/(m K) % Original: 0.15

par.d_wall = 0.03; %m

par.k_cargo = 120; %W/(m^2 K) -> Air / Aluminum
par.k_air_cooler = 120; %W/(m^2 K) -> Air / Aluminium (120)
par.k_air_wall = 15; %W/(m^2 K) -> Air / Wall
par.k_air_wall_i = 15; %W/(m^2 K) -> Wall inside

par.k_air_cooler_erzwungen = 150.101; %W/(m^2 K) -> Air / Aluminium
par.k_air_cooler_frei = 50; %W/(m^2 K) -> Air / Aluminium

% par.T_atm = 20+273.15; %Degree Celsius
% par.T_water = 20+273.15; %Degree Celsius

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
T_amb_mean = par.T_amb_mean;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
par.c = c;
% save('c.mat','c');



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

%% Iteration 1:

Parameters = [];

Parameters(1) = par.a1;
Parameters(2) = par.a2;
Parameters(3) = par.a3;
% Parameters(3) = par.a4;
Parameters(4) = par.a6;
Parameters(5) = par.a_frei_1;
Parameters(6) = par.a7;
Parameters(7) = par.b2;
Parameters(8) = par.a8;
Parameters(9) = par.a9;
Parameters(10) = par.a11;
Parameters(11) = par.a_frei_2;
Parameters(12) = par.xi;
Parameters(13) = par.b3;

for i = 1:numel(Parameters)
    if Parameters(i) < 0
        Parameters(i) = Parameters(i)*(-1);
    end
end

Parameters = Parameters';

Theta_1_float = 0.5;

f_1 = [];
f_1 = Parameters./Theta_1_float;

for i = 1:numel(Parameters)
    Theta_0(i) = Theta_1_float;
end
Theta_0 = Theta_0';

par.Parameters = Parameters;
par.Theta_0 = Theta_0;
par.f_1 = f_1;

% save('f_1_a4_new.mat','f_1');




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

x0_x1 = Data_18(1,3) + 273.15;
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
sol_comparison = ode15s(@(t,x) ICB_comparison(t,x,par), tspan, x0, options);

% Solution Values:
sim = (deval(sol,t_ref))';
sim(:,:) = sim(:,:)-273.15;
t_ref_h = t_ref./3600;

x_sim = (sol.y)';
t_sim = (sol.x)';
x_sim(:,1:5) = x_sim(:,1:5)-273.15;

% Comparison Values:
sim_comp = (deval(sol_comparison,t_ref))';
sim_comp(:,:) = sim_comp(:,:)-273.15;

x_sim_comp = (sol_comparison.y)';
t_sim_comp = (sol_comparison.x)';
x_sim_comp(:,1:5) = x_sim_comp(:,1:5)-273.15;



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

%% Save Parameters:

% First T_wall_inside:
T_wall_inside = sim(:,1); %In Degree Celsius

T_wall_inside_1(:,1) = T_wall_inside(1:3227,1);
T_wall_inside_2(:,1) = T_wall_inside(3228:(3227*2),1);
T_wall_inside_3(:,1) = T_wall_inside(((3227*2)+1):end,1);

% save('T_wall_inside.mat','T_wall_inside');
% save('T_wall_inside_1.mat','T_wall_inside_1');
% save('T_wall_inside_2.mat','T_wall_inside_2');
% save('T_wall_inside_3.mat','T_wall_inside_3');

% T_wall_outside:

T_wall_outside = sim(:,2); %In Degree Celsius

T_wall_outside_1(:,1) = T_wall_outside(1:3227,1);
T_wall_outside_2(:,1) = T_wall_outside(3228:(3227*2),1);
T_wall_outside_3(:,1) = T_wall_outside(((3227*2)+1):end,1);

%This is necessary for the Identification!
% save('T_wall_outside.mat','T_wall_outside');
% save('T_wall_outside_1.mat','T_wall_outside_1');
% save('T_wall_outside_2.mat','T_wall_outside_2');
% save('T_wall_outside_3.mat','T_wall_outside_3');


% save('Parameters.mat','par');



%% Plots




% Zusammenfassung der Simulation MIT DEVAL!:
% figure('Name','Simulation with Kelvin')
% yline(T_atm,'--')
% hold on
% plot(t_ref_h(:),sim(:,1),'Linewidth',1)
% plot(t_ref_h(:),sim(:,2),'Linewidth',1)
% plot(t_ref_h(:),sim(:,3),'Linewidth',1)
% plot(t_ref_h(:),sim(:,4),'Linewidth',1)
% plot(t_ref_h(:),sim(:,5),'Linewidth',1)
% hold off
% xlim([t_beginn,t_end/3600])
% ylim([round(min(sim(:,1:5),[],'all')-1,0),max(sim(:,1:5),[],'all')+5])
% % xticks([0, 0.5*3600, 3600, 1.5*3600, 2*3600, 2.5*3600, 3*3600, 3.5*3600,...
% %     4*3600, 4.5*3600, 5*3600, 5.5*3600, 6*3600, 6.5*3600, 7*3600, ...
% %     7.5*3600, 8*3600, 8.5*3600, 9*3600, 9.5*3600, 10*3600])
% % xticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5',...
% %     '5.5','6','6.5', '7', '7.5', '8', '8.5', '9', '9.5', '10'})
% xlabel('Time [h]')
% ylabel('Temperature [\circC]')
% legend('T_{atm}','T_{wall.inside}','T_{wall.outside}','T_{Cargo}','T_{HS}',...
%     'T_{ICB}')
% title('Simulation with Kelvin')


% Vergleich der ICB Temperatur:
figure('Name','Difference T_ICB')
plot(t_ref_h(:),sim(:,5),'Linewidth',1)
hold on
plot(t_ref_h(:),Data_18(:,3),'Linewidth',1)
hold off
xlim([t_beginn,t_end/3600])
ylim([round(min(sim(:,1:5),[],'all')-1,0),max(sim(:,1:5),[],'all')+5])
xlabel('Time [h]')
ylabel('Temperature [\circC]')
legend('T_{Simulation}','T_{Actual}')
title('Difference T_{ICB}')

% Vergleich der HS Temperatur:
figure('Name','Difference T_HS')
plot(t_ref_h(:),sim(:,4),'Linewidth',1)
hold on
plot(t_ref_h(:),Data_18(:,2),'Linewidth',1)
hold off
xlim([t_beginn,t_end/3600])
ylim([round(min(Data_18(:,2),[],'all')-1,0),max(sim(:,4),[],'all')+5])
xlabel('Time [h]')
ylabel('Temperature [\circC]')
legend('T_{Simulation}','T_{Actual}')
title('Difference T_{HS}')

% Vergleich der Wandinnen Temperatur:
% figure('Name','Difference T_Wall.Inside')
% plot(t_ref_h(:),sim(:,1),'Linewidth',1)
% hold on
% plot(t_ref_h(:),Data_18(:,4),'Linewidth',1)
% hold off
% xlim([t_beginn,t_end/3600])
% ylim([round(min(Data_18(:,4),[],'all')-1,0),max(sim(:,1),[],'all')+5])
% xlabel('Time [h]')
% ylabel('Temperature [\circC]')
% legend('T_{Simulation}','T_{Actual}')
% title('Difference T_{Wall.Inside}')



% % Zusammenfassung der gesamten Simulation OHNE DEVAL!
% figure('Name','Simulation with Kelvin')
% yline(T_atm,'--')
% hold on
% plot(t_sim(:),x_sim(:,1),'Linewidth',1)
% plot(t_sim(:),x_sim(:,2),'Linewidth',1)
% plot(t_sim(:),x_sim(:,3),'Linewidth',1)
% plot(t_sim(:),x_sim(:,4),'Linewidth',1)
% plot(t_sim(:),x_sim(:,5),'Linewidth',1)
% hold off
% xlim([t_beginn,t_end])
% ylim([round(min(x_sim(:,1:5),[],'all')-1,0),max(x_sim(:,1:5),[],'all')+5])
% xticks([0, 0.5*3600, 3600, 1.5*3600, 2*3600, 2.5*3600, 3*3600, 3.5*3600,...
%     4*3600, 4.5*3600, 5*3600, 5.5*3600, 6*3600, 6.5*3600, 7*3600, ...
%     7.5*3600, 8*3600, 8.5*3600, 9*3600, 9.5*3600, 10*3600])
% xticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5',...
%     '5.5','6','6.5', '7', '7.5', '8', '8.5', '9', '9.5', '10'})
% xlabel('Time [h]')
% ylabel('Temperature [\circC]')
% legend('T_{atm}','T_{wall.inside}','T_{wall.outside}','T_{Cargo}','T_{HS}',...
%     'T_{ICB}')
% title('Simulation with Kelvin')

%% Plot Comparison to evaluate:

% First an error Calculation is necessary:
Error_HS = sim(:,4)-sim_comp(:,4);
Error_ICB = sim(:,5)-sim_comp(:,5);


% Vergleich der HS Temperatur:

figure('Name','Errorfunction T_HS')
plot(t_ref_h(:),Error_HS(:),'Linewidth',1)
xlim([t_beginn,t_end/3600])
% ylim([round(min(sim_comp(:,4),[],'all')-1,0),max(sim(:,4),[],'all')+5])
xlabel('Time [h]')
ylabel('Temperature [\circC]')
% legend('T_{Simulation}','T_{Comparison}')
title('Error T_{HS}')

% figure('Name','Comparison T_HS')
% plot(t_ref_h(:),sim(:,4),'Linewidth',1)
% hold on
% plot(t_ref_h(:),sim_comp(:,4),'Linewidth',1)
% hold off
% xlim([t_beginn,t_end/3600])
% ylim([round(min(sim_comp(:,4),[],'all')-1,0),max(sim(:,4),[],'all')+5])
% xlabel('Time [h]')
% ylabel('Temperature [\circC]')
% legend('T_{Simulation}','T_{Comparison}')
% title('Comparison T_{HS}')

% Vergleich der ICB Temperatur:


figure('Name','Errorfunction T_ICB')
plot(t_ref_h(:),Error_ICB(:),'Linewidth',1)
xlim([t_beginn,t_end/3600])
% ylim([round(min(sim_comp(:,4),[],'all')-1,0),max(sim(:,4),[],'all')+5])
xlabel('Time [h]')
ylabel('Temperature [\circC]')
% legend('T_{Simulation}','T_{Comparison}')
title('Error T_{ICB}')


% figure('Name','Comparison T_ICB')
% plot(t_ref_h(:),sim(:,5),'Linewidth',1)
% hold on
% plot(t_ref_h(:),sim_comp(:,5),'Linewidth',1)
% hold off
% xlim([t_beginn,t_end/3600])
% ylim([round(min(sim(:,1:5),[],'all')-1,0),max(sim(:,1:5),[],'all')+5])
% xlabel('Time [h]')
% ylabel('Temperature [\circC]')
% legend('T_{Simulation}','T_{Comparison}')
% title('Comparison T_{ICB}')



%% Differential Equation (Function):

function dxdt_comparison = ICB_comparison(t,x,par)


%unpack all Parameters:
Parameters = par.Parameters;
Theta_0 = par.Theta_0;
f_1 = par.f_1;

th = Theta_0(1);

% T_atm = par.T_atm;
% T_water = par.T_water;

s_cargo = par.s_cargo;
c = par.c;


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
dxdt_comparison(1) = f_1(1)*th*(x(2)-x(1)) - f_1(2)*th*(x(1)-x(5));
dxdt_comparison(2) = f_1(3)*th*(T_amb-x(2)) - c(2)*(x(2)-x(1));
dxdt_comparison(3) = 0;
dxdt_comparison(4) = (f_1(4)*s_Fan+(1-s_Fan)*f_1(5))*th*(x(5)-x(4)) - ... 
    f_1(6)*th*I*x(4) + f_1(7)*th*(I^2) + f_1(8)*th*(T_water-x(4));
dxdt_comparison(5) = f_1(9)*th*(x(1)-x(5)) - ...
    (f_1(10)*s_Fan+(1-s_Fan)*f_1(11))*th*(x(5)-x(4)) + ...
    f_1(12)*th*(T_amb-x(5))*s_Door + f_1(13)*th*s_Fan;



dxdt_comparison = dxdt_comparison';

end


function dxdt = ICB(t,x,par)


%unpack all Parameters:
a1 = par.a1;
a2 = par.a2;
a3 = par.a3;
a4 = par.a4;
a5 = par.a5;
a6 = par.a6; % Includes the forced convection
a7 = par.a7;
a8 = par.a8;
a9 = par.a9;
a10 = par.a10;
a11 = par.a11; % Includes the forced convection

a_frei_1 = par.a_frei_1; % Free convection parameter 
a_frei_2 = par.a_frei_2; % Free convection parameter

% b1 = par.b1;
b2 = par.b2;
b3 = par.b3;


% T_atm = par.T_atm;
% T_water = par.T_water;

s_cargo = par.s_cargo;

c = par.c;


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

dxdt(1) = a1*(x(2)-x(1)) - a2*(x(1)-x(5));
dxdt(2) = a3*(T_amb-x(2)) - a4*(x(2)-x(1));
dxdt(3) = 0;
dxdt(4) = (a6*s_Fan + (1-s_Fan)*a_frei_1)*(x(5)-x(4)) - a7*I*x(4) + ...
    b2*(I^2) + a8*(T_water-x(4));
dxdt(5) = a9*(x(1)-x(5)) - ...
    (a11*s_Fan + (1-s_Fan)*a_frei_2)*(x(5)-x(4)) + ...
    xi*(T_amb-x(5))*s_Door + b3*s_Fan;

dxdt = dxdt';

end



