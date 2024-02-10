%% BSC-Thesis

%Author: Yannick Ramic
%Date: 27.07.2022
%Mtr-Nr.: 11771174


clear all;
clc;

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

par.T_atm = 20+273.15; %Degree Celsius ... Original: 20
par.T_water = 15+273.15; %Degree Celsius ... Original: 20

%Temperature of the wall outside can be described as a constant and has the
%same value than the ambient temperature
% par.T_wall_o = 20; 

par.P_fan = 0.5; %W
par.c_fan = 0.5; %W

par.I = 3; %A (Amper) -> 4.7, Used Value! 7.5302
par.I_open = 0; % 3, 6

%Peltier Element:

par.U_max = 20; %V
par.I_max = 8.5; %A
par.delta_T_max = 71; %K or Degree Celsius

%First all necessary Parameters for the constats alpha, R & K need to get
%unpacked:

U_max = par.U_max;
I_max = par.I_max;
delta_T_max = par.delta_T_max;
T_water = par.T_water;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_hot = T_water + 273.15;
% par.T_hot = T_hot;
% 
% alpha_PE = U_max / T_hot;
% R_PE = ((T_hot - delta_T_max) * U_max) / (T_hot * I_max);
% K_PE = ((T_hot - delta_T_max) * U_max * I_max) / (2 * T_hot * delta_T_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha_PE = U_max / T_water;
R_PE = ((T_water - delta_T_max) * U_max) / (T_water * I_max);
K_PE = ((T_water - delta_T_max) * U_max * I_max) / (2 * T_water * delta_T_max);

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

T_atm = par.T_atm;
T_water = par.T_water;

% T_wall_o = par.T_wall_o; 

P_fan = par.P_fan;
c_fan = par.c_fan;


I = par.I;





%% System Coefficients:

par.a1 = (A_wall*lambda_wall)/(m_wall*c_wall*d_wall);
par.a2 = (-A_wall*k_air_wall)/(m_wall*c_wall);
par.a3 = (k_air_wall*A_wall)/(m_wall*c_wall);
par.a4 = (-A_wall*lambda_wall)/(m_wall*c_wall*d_wall);
par.a5 = (-A_cargo*k_cargo)/(m_cargo*c_cargo);
par.a6 = (2*k_air_cooler*A_hs)/(m_hs*c_hs);
par.a7 = -2*alpha_PE/(m_hs*c_hs); %!!! Linearisierung
par.a8 = 2*K_PE/(m_hs*c_hs);
% par.9 siehe unten
par.a10 = (A_cargo*k_cargo)/(m_air*c_air);
par.a9 = (A_wall*k_air_wall)/(m_air*c_air);
par.a11 = (-2*k_air_cooler*A_hs)/(m_air*c_air);

% par.b1 = ...
par.b2 = 2*R_PE/(2*m_hs*c_hs);
par.b3 = 2*c_fan/(m_air*c_air);

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



%% Find x1.stat:

% Assumptions to find a x1.stat:
T_ICB_stat = 16 + 273.15; % 18
I_stat = 3;

%Formulate Equations:
syms T_wi_stat T_wo_stat T_cargo_stat T_HS_stat xi_stat
eqn_1 = a1*(T_wo_stat - T_wi_stat) + a2*(T_wi_stat - T_ICB_stat) == 0;
eqn_2 = a3*(T_atm - T_wo_stat) + a4*(T_wo_stat - T_wi_stat) == 0;
eqn_3 = a5*(T_cargo_stat - T_ICB_stat) == 0;
eqn_4 = a6*(T_ICB_stat - T_HS_stat) + a7*T_HS_stat*I_stat + ... 
    b2*(I_stat^2) + a8*(T_water-T_HS_stat) == 0;
eqn_5 = xi_stat*abs(T_atm - T_ICB_stat) + a9*(T_wi_stat - T_ICB_stat) + ...
    a10*0*(T_cargo_stat - T_ICB_stat) + a11*(T_ICB_stat - T_HS_stat) == 0;

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

par.c_D = 1; %Discharge Coefficient -> LaFayer Paper (0.85- 1!)

par.c1_v = (par.g * par.H) / (8*par.L);
par.c2_v = - 9 / (8*par.L*(par.c_D ^2));
par.c3_v = (2*par.W_AP) / (par.L*par.W);
par.c4_v = 1 / (par.m_air*par.c_air);
par.c5_v = par.c3_v / par.c4_v;

% x1.stat:
x1_stat_cal = Solution_xi_stat(5,1) / (2*par.c3_v);
par.x1_stat = x1_stat_cal; %See LaFayer Paper (x1.stat = 0.1)
% par.x1_stat = 0.15;
x1_stat = par.x1_stat;

par.xi = 2 * par.c3_v * par.x1_stat;
xi = par.xi;

par.s_door = 0; %Tür geschlossen
% par.a9 = -par.xi;
% par.e3 = -par.a9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Linearization of the system while Door closed:

% Assumptions to find a x1.stat:
T_ICB_lin = 16 + 273.15;
T_amb_lin = 20 + 273.15;
T_water_lin = 15 + 273.15;

%Formulate Equations:
syms T_wi_lin T_wo_lin T_cargo_lin T_HS_lin I_lin
eqn_1_lin = a1*(T_wo_lin - T_wi_lin) + a2*(T_wi_lin - T_ICB_lin) == 0;
eqn_2_lin = a3*(T_amb_lin - T_wo_lin) + a4*(T_wo_lin - T_wi_lin) == 0;
eqn_3_lin = a5*(T_cargo_lin - T_ICB_lin) == 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eqn_4_lin = a6*(T_ICB_lin - T_HS_lin) + a7*T_HS_lin*I_lin + ... 
%     b2*(I_lin^2) + a8*(T_water_lin - T_HS_lin) == 0;
% eqn_5_lin = xi*abs(T_amb_lin - T_ICB_lin) + a9*(T_wi_lin - T_ICB_lin) + ...
%     a10*(T_cargo_lin - T_ICB_lin) + a11*(T_ICB_lin - T_HS_lin) == 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eqn_5_lin = a9*(T_wi_lin - T_ICB_lin) + a10*(T_cargo_lin - T_ICB_lin) ...
    + a11*(T_ICB_lin - T_HS_lin) + b3 == 0;


% Simplify the system of Equations into a Matrix and a Vektor:

% [Matrix_A_lin, Vector_b_lin] = equationsToMatrix([eqn_1_lin, eqn_2_lin, ...
%     eqn_3_lin,eqn_4_lin, eqn_5_lin], [T_wi_lin, T_wo_lin, T_cargo_lin, ...
%     T_HS_lin, I_lin]);
[Matrix_A_lin, Vector_b_lin] = equationsToMatrix([eqn_1_lin, eqn_2_lin, ...
    eqn_3_lin, eqn_5_lin], [T_wi_lin, T_wo_lin, T_cargo_lin, ...
    T_HS_lin]);

% Solve the system of Equations:
Solution_lin = linsolve(Matrix_A_lin,Vector_b_lin);

Solution_lin_stat = [];
%Spalte 1 -> T_wi_lin
%Spalte 2 -> T_wo_lin
%Spalte 3 -> T_cargo_lin
%Spalte 4 -> T_HS_lin

for i = 1:4
    Solution_lin_stat(i,1) = Solution_lin(i,1)-273.15;  
end

Solution_lin_K = Solution_lin_stat + 273.15;

eqn_4_lin = a6*(T_ICB_lin - Solution_lin_K(4,1)) + ...
    a7*Solution_lin_K(4,1)*I_lin + b2*(I_lin^2) + ...
    a8*(T_water_lin - Solution_lin_K(4,1)) == 0;

% 2 Solutions for I_lin due to the fact that the equation is quadratic
I_lin_stat = solve(eqn_4_lin,I_lin);
I_lin_sol = [];
I_lin_sol(1:2,1) = I_lin_stat(1:2,1);
I_lin_closed = I_lin_sol(1,1);

lin_K = Solution_lin_K;
lin_K(5,1) = T_ICB_lin;

par.lin_K = lin_K;

%% FIND X1.STAT:

% Assumptions to find a x1.stat:
T_ICB_AP = 14 + 273.15;
T_amb_AP = 20 + 273.15;
T_water_AP = 15 + 273.15;

%Formulate Equations:
syms T_wi_AP T_wo_AP T_cargo_AP T_HS_AP I_AP
eqn_1_AP = a1*(T_wo_AP - T_wi_AP) + a2*(T_wi_AP - T_ICB_AP) == 0;
eqn_2_AP = a3*(T_amb_AP - T_wo_AP) + a4*(T_wo_AP - T_wi_AP) == 0;
eqn_3_AP = a5*(T_cargo_AP - T_ICB_AP) == 0;

eqn_5_AP = a9*(T_wi_AP - T_ICB_AP) + a10*(T_cargo_AP - T_ICB_AP) + ...
    a11*(T_ICB_AP - T_HS_AP) + xi*(T_amb_AP - T_ICB_AP) == 0;

[Matrix_A_AP, Vector_b_AP] = equationsToMatrix([eqn_1_AP, eqn_2_AP, ...
    eqn_3_AP, eqn_5_AP], [T_wi_AP, T_wo_AP, T_cargo_AP, T_HS_AP]);

% Solve the system of Equations:
Solution_AP = linsolve(Matrix_A_AP,Vector_b_AP);

Solution_AP_stat = [];
%Spalte 1 -> T_wi_AP
%Spalte 2 -> T_wo_AP
%Spalte 3 -> T_cargo_AP
%Spalte 4 -> T_HS_AP

for i = 1:4
    Solution_AP_stat(i,1) = Solution_AP(i,1)-273.15;  
end

Solution_AP_K = Solution_AP_stat + 273.15;

eqn_4_AP = a6*(T_ICB_AP - Solution_AP_K(4,1)) + ...
    a7*Solution_AP_K(4,1)*I_AP + b2*(I_AP^2) + ...
    a8*(T_water_AP - Solution_AP_K(4,1)) == 0;

% 2 Solutions for I_lin due to the fact that the equation is quadratic
I_AP_stat = solve(eqn_4_AP,I_AP);
I_AP_sol = [];
I_AP_sol(1:2,1) = I_AP_stat(1:2,1);
I_open = I_AP_sol(1,1);

I_open = 0; % ADDED!!!


lin_K_open = Solution_AP_K;
lin_K_open(5,1) = T_ICB_AP;

par.lin_K_open = lin_K_open;

%% Initial Condition Problem + Numerical Coefficients:

x0_x1 = 20 + 273.15;
x0_x2 = 20 + 273.15;
x0_x3 = 20 + 273.15;
x0_x4 = 20 + 273.15;
x0_x5 = 20 + 273.15;
x0_x6 = 0;
x0 = [x0_x1; x0_x2; x0_x3; x0_x4; x0_x5; x0_x6];

%Time:

%All Time Variables are in seconds
t_step = 3*3600; 
t_open = 6*3600;
t_close = (t_open + 150);
t_beginn = 0;
t_end = 10*3600; 
tspan = [t_beginn,t_end];

par.t_step = t_step;
par.t_open = t_open;
par.t_close = t_close;

%% ODE Solver:

options = odeset('Events', @(t,x) Event(t,x,par));
sol = ode15s(@(t,x) icb(t,x,par), tspan, x0, options);

sol_1 = sol;

counter = 1;
I_1 = par.I;
I_2 = 3; % 7.5302
par.I = I_2;

I = [];
x = [];
t = [];

x0_step = sol.y(:,end);
sol = odextend(sol, @(t,x) icb(t,x,par), t_end, x0_step, options);

x0_open = sol.y(:,end);
sol = odextend(sol, @(t,x) open_door(t,x,par), t_end, x0_open, options);

x0_close = sol.y(:,end);
x0_close(6,1) = x0(6,1);
sol = odextend(sol, @(t,x) icb(t,x,par), t_end, x0_close, options);

while sol.x(counter) < t_end
    
    if sol.x(counter) <= t_step
        
        I(counter,1)=sol.x(counter);
        I(counter,2)= I_1;
        
        x(counter,:) = (sol.y(:,counter))';
        t(counter) = (sol.x(counter))';
        
    else
        
        I(counter,1)=sol.x(counter);
        I(counter,2) = I_2;
        
    end
   
    counter = counter + 1;
    
end


x_sim = (sol.y)';
t_sim = (sol.x)';
x_sim(:,1:5) = x_sim(:,1:5)-273.15;

%% Linearized Model for Door Closed:

% The shape of Vektor x as followed:
% Whereas the index 0 stands for a point of interest:
% 1) T_w_i - T_w_i_0
% 2) T_w_o - T_w_o_0
% 3) T_cargo - T_cargo_0
% 4) T_HS - T_HS_0
% 5) T_ICB - T_ICB_0

Zustand = [(a2-a1), a1, 0, 0, -a2;...
    -a4, (a4-a3), 0, 0, 0;...
    0, 0, a5, 0, -a5;...
    0, 0, 0, (a7*I_lin_closed-a6-a8), a6;...
    a9, 0, a10, -a11, (a11-a10-a9)];

% Shape of u (Input Values):
% 1) I - I_0
% 2) s_fan - s_fan_0 = 0
% 1) s_door - s_door_0 = 0

Input = [0, 0, 0;...
    0, 0, 0;...
    0, 0, 0;...
    (a7*Solution_lin_K(4,1)+2*b2*I_lin_closed), 0, 0;...
    0, b3, 0];

% Shape of z (Disturbances):
% 1) T_amb - T_amb_0 = 0 ... (25 - 25)
% 2) T_water - T_water_0 = 0 ... (15 - 15)

Disturbance = [0, 0;...
    a3, 0;...
    0, 0;...
    0, a8;...
    0, 0];


delta_u = [I_1 - I_lin_closed; 0; 0];
delta_z = [0; 0];

Disturbance_lin = Disturbance*delta_z;
Input_lin = Input*delta_u;

par.Disturbance_lin = Disturbance_lin;
par.Input_lin = Input_lin;
par.Zustand = Zustand;


%% Comparison of the non linear to the linearized model (Door Closed):

% sol_lin = ode15s(@(t,x) icb_lin(t,x,par), [t_beginn, t_step], x0(1:5,1));
% 
% eval = linspace(0, t_step);
% sol_lin_eval = deval(sol_lin, eval);
% 
% x_lin = (sol_lin_eval)';
% t_lin = eval';
% 
% % x_lin = (sol_lin.y)';
% % t_lin = (sol_lin.x)';
% % x_lin(:,1:5) = x_lin(:,1:5)-273.15;
% 
% t_diff = [];
% x_diff = [];

% for i = 1:numel(t_sim)
%     if t_sim(i) <= t_step
%         t_diff(i,1) = t_sim(i);
%         x_diff(i,1:5) = x_sim(i,1:5);
%     end
% end

% x_diff_eval = deval(sol, eval)';
% x_diff_eval(:,1:5) = x_diff_eval(:,1:5)-x_lin;

% figure('Name','Difference')
% plot(t_lin,x_diff_eval(:,1:5),'Linewidth',1)
% xlim([t_beginn,t_step])
% % ylim([round(min(x_diff_eval(:,1:5),[],'all'),0),max(x_diff_eval(:,1:5),[],'all')+0.25])
% xticks([0, 0.5*3600, 3600, 1.5*3600, 2*3600, 2.5*3600, 3*3600])
% xticklabels({'0','0.5','1','1.5','2','2.5','3'})
% xlabel('Time [h]')
% ylabel('Temperature \DeltaT [\circC]')
% legend('\DeltaT_{wall.inside}','\DeltaT_{wall.outside}',...
%     '\DeltaT_{Cargo}','\DeltaT_{HS}', '\DeltaT_{ICB}')
% title('Differences due to Linearization')
% 
% x_lin(:,1:5) = x_lin(:,1:5)-273.15;
% 
% figure('Name','Linearized Model')
% % yline(T_atm,'--')
% % hold on
% plot(t_lin(:),x_lin(:,1),'Linewidth',1)
% hold on
% plot(t_lin(:),x_lin(:,2),'Linewidth',1)
% plot(t_lin(:),x_lin(:,3),'Linewidth',1)
% plot(t_lin(:),x_lin(:,4),'Linewidth',1)
% plot(t_lin(:),x_lin(:,5),'Linewidth',1)
% % xline(t_step)
% % xline(t_open)
% % xline(t_close)
% hold off
% xlim([t_beginn,t_step])
% ylim([round(min(x_lin(:,1:5),[],'all')-1,0),max(x_lin(:,1:5),[],'all')+5])
% % xticks([0, 0.5*3600, 3600, 1.5*3600, 2*3600, 2.5*3600, 3*3600, 3.5*3600,...
% %     4*3600, 4.5*3600, 5*3600, 5.5*3600, 6*3600, 6.5*3600, 7*3600, ...
% %     7.5*3600, 8*3600, 8.5*3600, 9*3600, 9.5*3600, 10*3600])
% % xticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5',...
% %     '5.5','6','6.5', '7', '7.5', '8', '8.5', '9', '9.5', '10'})
% xlabel('Time [h]')
% ylabel('Temperature [\circC]')
% % legend('T_{atm}','T_{wall.inside}','T_{wall.outside}','T_{Cargo}','T_{HS}',...
% %     'T_{ICB}')
% legend('T_{wall.inside}','T_{wall.outside}','T_{Cargo}','T_{HS}', 'T_{ICB}')
% title('Linearized Model')

%% Linearized Model for Door Open:

% The shape of Vektor x as followed:
% Whereas the index 0 stands for a point of interest:
% 1) T_w_i - T_w_i_0
% 2) T_w_o - T_w_o_0
% 3) T_cargo - T_cargo_0
% 4) T_HS - T_HS_0
% 5) T_ICB - T_ICB_0

Zustand_open = [(a2-a1), a1, 0, 0, -a2;...
    -a4, (a4-a3), 0, 0, 0;...
    0, 0, a5, 0, -a5;...
    0, 0, 0, (a7*I_open-a6-a8), a6;...
    a9, 0, a10, -a11, (-xi+a11-a10-a9)];

% Shape of u (Input Values):
% 1) I - I_0
% 2) s_fan - s_fan_0 = 0
% 1) s_door - s_door_0 = 0

Input_open = [0, 0, 0;...
    0, 0, 0;...
    0, 0, 0;...
    (a7*Solution_AP_K(4,1)+2*b2*I_open), 0, 0;...
    0, b3, xi*(T_atm-T_ICB_AP)];

% Shape of z (Disturbances):
% 1) T_amb - T_amb_0 = 0 ... (25 - 25)
% 2) T_water - T_water_0 = 0 ... (15 - 15)

Disturbance_open = [0, 0;...
    a3, 0;...
    0, 0;...
    0, a8;...
    xi, 0];


delta_u_open = [0; 0; 0];
delta_z_open = [0; 0];

Disturbance_lin_open = Disturbance_open*delta_z_open;
Input_lin_open = Input_open*delta_u_open;

par.Disturbance_lin_open = Disturbance_lin_open;
par.Input_lin_open = Input_lin_open;
par.Zustand_open = Zustand_open;

sol_lin = ode15s(@(t,x) icb_lin(t,x,par), tspan, x0(1:5,1),options);

x0_step_lin = sol_lin.y(:,end);
sol_lin = odextend(sol_lin, @(t,x) icb_lin(t,x,par), t_end, x0_step_lin, options);

x0_open_lin = sol_lin.y(:,end);
sol_lin = odextend(sol_lin, @(t,x) icb_lin_open(t,x,par), t_end, x0_open_lin, options);

x0_close_lin = sol_lin.y(:,end);
sol_lin = odextend(sol_lin, @(t,x) icb_lin(t,x,par), t_end, x0_close_lin, options);

x_lin_sim = sol_lin.y' - 273.15;
t_lin_sim = sol_lin.x';


figure('Name','Linearized Model')
yline(T_atm,'--')
hold on
plot(t_lin_sim,x_lin_sim,'Linewidth',1)
xline(t_step)
xline(t_open)
xline(t_close)
hold off
xlim([t_beginn,t_end])
ylim([round(min(x_lin_sim(:,1:5),[],'all')-1,0),max(x_lin_sim(:,1:5),[],'all')+5])
xticks([0, 0.5*3600, 3600, 1.5*3600, 2*3600, 2.5*3600, 3*3600, 3.5*3600,...
    4*3600, 4.5*3600, 5*3600, 5.5*3600, 6*3600, 6.5*3600, 7*3600, ...
    7.5*3600, 8*3600, 8.5*3600, 9*3600, 9.5*3600, 10*3600])
xticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5',...
    '5.5','6','6.5', '7', '7.5', '8', '8.5', '9', '9.5', '10'})
xlabel('Time [h]')
ylabel('Temperature [\circC]')
% legend('T_{atm}','T_{wall.inside}','T_{wall.outside}','T_{Cargo}','T_{HS}',...
%     'T_{ICB}')
legend('T_{wall.inside}','T_{wall.outside}','T_{Cargo}','T_{HS}', 'T_{ICB}')
title('Linearized Model')

%% Heat Analysis during Door Opening

c5_v = par.c5_v;

i_open = find(x_sim(:,6) ~= 0);
delta_T_icb = x_sim(i_open,5);
delta_t_open = t_sim(i_open);
delta_t_open_0 = delta_t_open - delta_t_open(1);

Q_door = (-2 * c5_v * x1_stat) * delta_T_icb + (2 * c5_v * x1_stat * ...
    (T_atm-273.15));

%% Store relevant Data:

Q_door_max = max(Q_door);

%First we need to find the indices of the door openings:
i_open = find(x_sim(:,6) ~= 0);
velocity = t_sim(i_open);
%Now I want that the velocity time starts from zero:
velocity = velocity - velocity(1,1);
velocity(:,2) = x_sim(i_open,6);

% Store Elements:

par_d.Q_door_max = Q_door_max;
par_d.xi = xi;
par_d.x1_stat = x1_stat;

par_CD1_T18 = par_d;
% save('par_CD1_T18.mat','par_CD1_T18')

Data_CD1_T18(1,:) = delta_t_open_0(:); % 1. Zeile: Time für Q
Data_CD1_T18(2,:) = Q_door(:); % 2. Zeile: Q_punkt Tür
Data_CD1_T18(3,:) = velocity(:,1); % 3. Zeile: Time für die Geschw.
Data_CD1_T18(4,:) = velocity(:,2); % 3. Zeile: Time für die Geschw.

% save('Data_CD1_T18.mat','Data_CD1_T18')


%% Plots

% First I checked the Heat Flow during the door opening whether the values
% make any sense. With a peak smaller than 700W and a door opening under 10
% min.

figure('Name','Heat Flow')
plot(delta_t_open_0,Q_door,'-o','Linewidth',1)
xlim([0, t_close-t_open])
ylim([0,max(Q_door,[],'all')+100])
xlabel('Time [s]')
ylabel('Q [W]')
title('Heat Flow During Door Opening Q_{Door}')



figure('Name','Simulation with Kelvin')
% yline(T_atm,'--')
hold on
plot(t_sim(:),x_sim(:,1),'Linewidth',1)
plot(t_sim(:),x_sim(:,2),'Linewidth',1)
plot(t_sim(:),x_sim(:,3),'Linewidth',1)
plot(t_sim(:),x_sim(:,4),'Linewidth',1)
plot(t_sim(:),x_sim(:,5),'Linewidth',1)
xline(t_step)
xline(t_open)
xline(t_close)
hold off
xlim([t_beginn,t_end])
ylim([round(min(x_sim(:,1:5),[],'all')-1,0),max(x_sim(:,1:5),[],'all')+5])
xticks([0, 0.5*3600, 3600, 1.5*3600, 2*3600, 2.5*3600, 3*3600, 3.5*3600,...
    4*3600, 4.5*3600, 5*3600, 5.5*3600, 6*3600, 6.5*3600, 7*3600, ...
    7.5*3600, 8*3600, 8.5*3600, 9*3600, 9.5*3600, 10*3600])
xticklabels({'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5',...
    '5.5','6','6.5', '7', '7.5', '8', '8.5', '9', '9.5', '10'})
xlabel('Time [h]')
ylabel('Temperature [\circC]')
legend('T_{wall.inside}','T_{wall.outside}','T_{Cargo}','T_{HS}',...
    'T_{ICB}')
title('Simulation with Kelvin')

% %First we need to find the indices of the door openings:
% i_open = find(x_sim(:,6) ~= 0);
% velocity = t_sim(i_open);
% %Now I want that the velocity time starts from zero:
% velocity = velocity - velocity(1,1);
% velocity(:,2) = x_sim(i_open,6);

figure('Name','Open_Door_Velocity')
plot(velocity(:,1),velocity(:,2),'Linewidth',1)
hold on
yline(x1_stat,'b--','Linewidth',1)
hold off
xlim([velocity(1,1), velocity(end,1)])
ylim([0,max(velocity(:,2))+0.05])
legend('$\dot{V}/A_{ap}$','$x_{1.stat}$','Interpreter','latex')
xlabel('Time [s]')
ylabel('$\dot{V}/A_{ap}$ [m/s]','Interpreter','latex')
title('Linearised Stream During Door Opening')

%% Differential Equation (Function):

%This Section needs to be at the end of the file!

function dxdt = icb(t,x,par)


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


% b1 = par.b1;
b2 = par.b2;
b3 = par.b3;


T_atm = par.T_atm;
T_water = par.T_water;

I = par.I;

%Beschreibung des Vektors x:
% x(1) -> T_wall_inside
% x(2) -> T_wall_outside
% x(3) -> T_cargo
% x(4) -> T_hs
% x(5) -> T_ICB
% x(6) -> (V_punkt / A_ap) [m/s]

%All 6 Differential Equations
dxdt(1) = a1*(x(2)-x(1)) + a2*(x(1)-x(5));
dxdt(2) = a3*(T_atm-x(2)) + a4*(x(2)-x(1));
dxdt(3) = a5*(x(3)-x(5));
dxdt(4) = a6*(x(5)-x(4)) + a7*I*x(4) + b2*(I^2) + a8*(T_water-x(4));
dxdt(5) = a9*(x(1)-x(5)) + a10*(x(3)-x(5)) + a11*(x(5)-x(4)) + b3*1;
dxdt(6) = 0;

dxdt = dxdt';

end


%This Function is an Event that stops the Iteration at some time which is
%set as a parameter.

function [value,isterminal,direction] = Event(t,~,par)

t_step = par.t_step;
t_open = par.t_open;
t_close = par.t_close;

value = [t-t_step; t-t_open; t-t_close]; 

isterminal = [1;1;1]; %terminate integration when y==0
direction = [0;0;0]; % 0-stellen bei: 1 steigender, -1 fallender, 0 jeder Funktionsstelle

end

%Function to solve the the ODE System for the open door:


function dxdt_o = open_door(t,x,par)


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


% b1 = par.b1;
b2 = par.b2;
b3 = par.b3;


T_atm = par.T_atm;
T_water = par.T_water;

I = par.I_open;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1_v = par.c1_v;
c2_v = par.c2_v;
c3_v = par.c3_v;
c4_v = par.c4_v;
x1_stat = par.x1_stat;
xi = par.xi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_water = par.T_water;

%Beschreibung des Vektors x:
% x(1) -> T_wall_inside
% x(2) -> T_wall_outside
% x(3) -> T_cargo
% x(4) -> T_hs
% x(5) -> T_ICB
% x(6) -> (V_punkt / A_ap) [m/s]


%All 6 Differential Equations
dxdt_o(1) = a1*(x(2)-x(1)) + a2*(x(1)-x(5));
dxdt_o(2) = a3*(T_atm-x(2)) + a4*(x(2)-x(1));
dxdt_o(3) = a5*(x(3)-x(5));
dxdt_o(4) = a6*(x(5)-x(4)) + a7*I*x(4) + b2*(I^2) + a8*(T_water-x(4));
dxdt_o(5) = a9*(x(1)-x(5)) + a10*(x(3)-x(5)) + a11*(x(5)-x(4)) + ...
    xi*(T_atm-x(5));
dxdt_o(6) = (4*c2_v * x1_stat) * x(6) + (-c1_v / T_atm)* x(5) + ...
    (-4 * c2_v * x1_stat^2 + c1_v);

dxdt_o = dxdt_o';

end

%Linearized Model for closed Door:

function dxdt_lin = icb_lin(t,x,par)

%unpack all Parameters:
Zustand = par.Zustand;
Disturbance_lin = par.Disturbance_lin;
Input_lin = par.Input_lin;
lin_K = par.lin_K;

% x(1) -> T_wall_inside
% x(2) -> T_wall_outside
% x(3) -> T_cargo
% x(4) -> T_hs
% x(5) -> T_ICB

% dxdt_lin = [];

%All 5 Differential Equations
Vektor_x = [x(1)-lin_K(1,1); x(2)-lin_K(2,1); x(3)-lin_K(3,1); ...
    x(4)-lin_K(4,1); x(5)-lin_K(5,1)];


for i = 1:numel(Input_lin)
    
    dxdt_lin(i) = Zustand(i,:)*Vektor_x + Input_lin(i) + Disturbance_lin(i);    
    
end

dxdt_lin = dxdt_lin';

end


function dxdt_lin_open = icb_lin_open(t,x,par)

%unpack all Parameters:
Zustand_open = par.Zustand_open;
Disturbance_lin_open = par.Disturbance_lin_open;
Input_lin_open = par.Input_lin_open;
lin_K_open = par.lin_K_open;

% x(1) -> T_wall_inside
% x(2) -> T_wall_outside
% x(3) -> T_cargo
% x(4) -> T_hs
% x(5) -> T_ICB

% dxdt_lin = [];

%All 5 Differential Equations
Vektor_x_open = [x(1)-lin_K_open(1,1); x(2)-lin_K_open(2,1);...
    x(3)-lin_K_open(3,1); x(4)-lin_K_open(4,1); x(5)-lin_K_open(5,1)];


for i = 1:numel(Input_lin_open)
    
    dxdt_lin_open(i) = Zustand_open(i,:)*Vektor_x_open +...
        Input_lin_open(i) + Disturbance_lin_open(i);    
    
end

dxdt_lin_open = dxdt_lin_open';

end


