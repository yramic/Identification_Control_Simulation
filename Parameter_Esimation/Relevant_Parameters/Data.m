%% DATEN FÜR DIE IDENTIFICATION IN EINEM FILE
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


close all; clear all; clc;

global c
load('c.mat')

load('Experiment_18.mat');
load('PARAMETER.mat');
load('T_wall_inside.mat');
load('T_wall_inside_1.mat');
load('T_wall_inside_2.mat');
load('T_wall_inside_3.mat');
load('T_wall_outside.mat');
load('T_wall_outside_1.mat');
load('T_wall_outside_2.mat');
load('T_wall_outside_3.mat');
load('Data_Set_1.mat');
load('Data_Set_2.mat');
load('Data_Set_3.mat');


% Aufbau Experiment_18:
% 1.Spalte: Time t
% 2.Spalte: T_ICB_actual WRONG!!!
% 3.Spalte: T_HS_actual
% 4.Spalte: T_ICB_actual
% 5.Spalte: I_actual
% 6.Spalte: s_Door (Door Status)
% 7.Spalte: s_Fan (Fan Status)
% 8.Spalte: I_ref_actual


% Parmeterfestlegung
% ------------------------------------------------------------------------

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

xi = par.xi;

% "unbekannte", konstante Systemparameter
theta = [a1;...
         a2;...
         a3;...
         a4;...
         a6;...
         a7;...
         b2;...
         a8;...
         a9;...
         a11;...
         xi;...
         b3];
     
 
% Gleichgewichtszustand definieren
% ------------------------------------------------------------------------
% Dieser Zustand wird im Folgenden als Ausgangsbasis für
% Stell-/Störgrößenveränderungen genutzt
x_stat = [T_wall_inside(1,1)+273.15, T_wall_outside(1,1)+273.15, ...
    M(1,3)+273.15, M(1,4)+273.15];
u_stat = [0.1; theta(3)*sqrt(x_stat(3))];
% d_stat = [];
% z_stat = [u_stat(1);0];

% Eingangs- und Störgrößen für einen Simulationsdurchgang festlegen
% ------------------------------------------------------------------------
t_vec = M(:,1);

% Eingangsgrößen
u(1,:) = M(:,8);
u(2,:) = M(:,6);
u(3,:) = M(:,7);

% messbare, nicht beeinflussbsare Eingangsgrößen (hier nicht vorhanden)
% d = [];

% Störgrößen
% z = zeros(2,size(t_vec,2));
% z(1,:) = z_stat(1,1);
% z(2,:) = z_stat(2,1);

% Startzustand für die Simulation soll der stationäre Zustand (_stat) sein
x_init = x_stat;

% Simulation durchführen und Ausgangsgrößen festlegen
% ------------------------------------------------------------------------
options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',2.5);
[t,x_num] = ode15s(@(t,x) Non_Linear_ICB(t,x,c,theta,u,t_vec),...
    [M(1,1) M(end,1)], x_init, options);

T_W_i = x_num(:,1)-273.15;
T_W_o = x_num(:,2)-273.15;
T_HS = x_num(:,3)-273.15;
T_ICB = x_num(:,4)-273.15;

% Simulation durchführen und Ausgangsgrößen festlegen
% ------------------------------------------------------------------------
figure
hold on
grid on, grid minor
plot(t,T_W_i);
plot(t,T_W_o);
plot(t,T_HS);
plot(t,T_ICB);
xlabel('Zeit in s');
ylabel('Temperatur in \circC');
legend('T_{Wall.Inside}','H_{Wall.Outside}','H_{HS}','H_{ICB}');

% Zusammenfassen der Eingangs- und Ausgangsgrößen in eine Matrix
% ------------------------------------------------------------------------
% Die Ausgangsgrößen (Füllstandshöhen) werden zur besseren Abbildung der
% Realität noch verrauscht.

% um eine zeitlich äquidistante Abtastung zu erhalten, wird hier noch
% zwischen den aus der Simulation erhaltenen Resultaten linear
% interpoliert.
% x_measure_num = interp1(t,x_num,t_vec); 
% 
% data = zeros(8,size(t_vec,2));
% data(1,:) = t_vec;
% data([2,3],:) = u;
% data([4,5],:) = z;
% data([6:8],:) = awgn(x_measure_num',60,'measured');

% optionales Abspeichern der erzeugten Messdaten in Messdatenordner
% save('01_measurement_data/Messung_3','data')


%% Global C:

T_atm = par.T_atm;
T_water = par.T_water;
c(1,1) = T_atm;
c(2,1) = T_water;

% save('c.mat','c')

%% Data Configuration 1:

% load('Messungen_1.mat');
% load('Messungen_1_2.mat');
% data = Messungen_1;
% 
% save('Messungen_1.mat','data')
% save('Messungen_1_2.mat','data')

%% Data Configuration 2:
% 
% load('Messungen_2.mat');
% load('Messungen_2_2.mat');
% data = Messungen_2;
% % 
% save('Messungen_2.mat','data')
% save('Messungen_2_2.mat','data')

%% Data Configuration 3:

% load('Messungen_3.mat');
% load('Messungen_3_2.mat');
% data = Messungen_3;
% 
% save('Messungen_3.mat','data')
% save('Messungen_3_2.mat','data')
