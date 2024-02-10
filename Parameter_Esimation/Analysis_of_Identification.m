%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Identification Process%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;



%% Set Up:

global c
global f_1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For more than one Optimization:
% global f_2
% load('f_2_a4.mat');


load('Experiment_18.mat');
load('Parameters.mat');
load('f_1_a4.mat');

load('c.mat');
% c(2) = 3.75e-4;
% c(2) = 8.346e-17;
addpath('Messungen_18');


data_1 = load('Messungen_1.mat');
data_1 = data_1.data;
data_2 = load('Messungen_2.mat');
data_2 = data_2.data;
data_3 = load('Messungen_3.mat');
data_3 = data_3.data;

data_tot = [data_1, data_2, data_3];

model_ident = load('model_a4_open.mat');
model = model_ident.model;
% save('model_a4_open_copy.mat','model')



% % % Aufbau Messungen_1
% % 
% % % 1.Zeile:          Zeit t
% % % 2.Zeile:          u(1) = I 
% % % 3.Zeile:          u(2) = s_Door
% % % 4.Zeile:          u(3) = s_Fan
% % % 5.Zeile:          x(1) = T_Wall_Inside
% % % 6.Zeile:          x(2) = T_Wall_Outside  
% % % 7.Zeile:          x(3) = T_HS
% % % 8.Zeile:          x(4) = T_ICB 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.Zeile:          Zeit t
% 2.Zeile:          u(1) = I 
% 3.Zeile:          u(2) = s_Door
% 4.Zeile:          u(3) = s_Fan
% 5.Zeile:          u(4) = T_water
% 6.Zeile:          u(5) = T_amb
% 7.Zeile:          x(1) = T_Wall_Inside
% 8.Zeile:          x(2) = T_Wall_Outside  
% 9.Zeile:          x(3) = T_HS
% 10.Zeile:         x(4) = T_ICB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% NICHTLINEARE GREY-BOX-IDENTIFIKATION
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% Diese Script führt eine Greybox-Identifikation des Kühlraumes durch

% Sampling time der Messdaten
Ts = 2.5;

%% VORBEREITUNG DER DATEN
% ************************************************************************
% Bezeichnung der Messungen

t = data_tot(1,1:end);
u = data_tot(2:6,1:end);
y = data_tot(9:10,1:end);
measurement = iddata(data_tot([9, 10],1:end)',u',Ts);

measurement.ExperimentName = 'Measurment';
measurement.InputName = {'u_1','u_2','u_3','u_4','u_5'};
measurement.InputUnit = {'A','1','1','K','K'};
measurement.OutputName = {'T_{HS}', 'T_{ICB}'};
measurement.OutputUnit = {'K','K'};
    
    

%% IDENTIFIKATION
% ************************************************************************
% Initialisierung mittels der bis dato angenommenen Systemparameter
% ------------------------------------------------------------------------

load('Data_Set_1.mat');
t_vec(1,:) = Data_Set_1(:,1);

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

% theta_1_0 = a1;
% theta_2_0 = a2;   
% theta_3_0 = a3;   
% theta_4_0 = a4;  
% theta_5_0 = a6;
% theta_6_0 = a7;   
% theta_7_0 = b2;
% theta_8_0 = a8;   
% theta_9_0 = a9;
% theta_10_0 = a11;
% theta_11_0 = xi;
% theta_12_0 = b3;

%All Parameter start with 0.5:

theta_0 = 0.5;

theta_1_0 = theta_0;
theta_2_0 = theta_0;   
theta_3_0 = theta_0;   
theta_4_0 = theta_0;  
theta_5_0 = theta_0;
theta_6_0 = theta_0;   
theta_7_0 = theta_0;
theta_8_0 = theta_0;   
theta_9_0 = theta_0;
theta_10_0 = theta_0;
theta_11_0 = theta_0;
theta_12_0 = theta_0;
theta_13_0 = theta_0;
% theta_14_0 = theta_0;

% Parameter festlegen, die im nichtlinearen Modell auftreten
% ------------------------------------------------------------------------
parameters = {theta_1_0; theta_2_0; theta_3_0; theta_4_0; theta_5_0;...
    theta_6_0; theta_7_0; theta_8_0; theta_9_0; theta_10_0; theta_11_0; ...
    theta_12_0; theta_13_0}; % BEFORE: theta_13_0; theta_14_0
Order = [2, 5, 4];      % Anzahl: Ausgang, Eingang, Zustände

InitialStates = [data_tot(7,1); data_tot(8,1); data_tot(9,1); data_tot(10,1)];

opt_c = compareOptions('InitialCondition','e'); % 'e' for estimated

figure
est = idnlgrey('Non_Linear_ICB_IDENT', Order, parameters, InitialStates);
est.SimulationOptions.Solver = 'ode15s';
compare(measurement,est,opt_c); 



% System zur Identifikation erzeugen
% ------------------------------------------------------------------------
% Es wird mittels eines kontinuierlichen Modells gearbeitet
sys = idnlgrey('Non_Linear_ICB_IDENT', Order, parameters, InitialStates);
sys.SimulationOptions.Solver = 'ode15s';


%%

Model_Data = model;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vergleich zwischen optimiertem Modell und Messung
% ------------------------------------------------------------------------
% if length(size(measurement_ident)) == 4
%     
%     num_exp_ident = size(measurement_ident); %Added
%     num_exp_ident = num_exp_ident(4); %Added
%     
%     for zaehler3 = 1:num_exp_ident
%         figure
%         compare(getexp(measurement_ident,zaehler3),model, opt_c);
%     end
% else
%     compare(measurement_ident,model, opt_c); 
% end

figure
compare(measurement,Model_Data, opt_c); 


clear Model_Data



%%
% Geschaetzte Parameter ausgeben
% ------------------------------------------------------------------------
disp('GESCHAETZTE PARAMETER:')
disp('--------------------------------------------------')
theta_1_mean = model.Report.Parameters.ParVector(1);
theta_1_std = sqrt(model.Report.Parameters.FreeParCovariance(1,1));
disp(['theta_1 = ', num2str(theta_1_mean), ' +/- ', num2str(theta_1_std), ' 1/s'])

theta_2_mean = model.Report.Parameters.ParVector(2);
theta_2_std = sqrt(model.Report.Parameters.FreeParCovariance(2,2));
disp(['theta_2 = ', num2str(theta_2_mean), ' +/- ', num2str(theta_2_std), ' 1/s'])

theta_3_mean = model.Report.Parameters.ParVector(3);
theta_3_std = sqrt(model.Report.Parameters.FreeParCovariance(3,3));
disp(['theta_3 = ', num2str(theta_3_mean), ' +/- ', num2str(theta_3_std), '1/s'])

theta_4_mean = model.Report.Parameters.ParVector(4);
theta_4_std = sqrt(model.Report.Parameters.FreeParCovariance(4,4));
disp(['theta_4 = ', num2str(theta_4_mean), ' +/- ', num2str(theta_4_std), ' 1/s'])

theta_5_mean = model.Report.Parameters.ParVector(5);
theta_5_std = sqrt(model.Report.Parameters.FreeParCovariance(5,5));
disp(['theta_5 = ', num2str(theta_5_mean), ' +/- ', num2str(theta_5_std), ' 1/s'])

theta_6_mean = model.Report.Parameters.ParVector(6);
theta_6_std = sqrt(model.Report.Parameters.FreeParCovariance(6,6));
disp(['theta_6 = ', num2str(theta_6_mean), ' +/- ', num2str(theta_6_std), ' V/J'])

theta_7_mean = model.Report.Parameters.ParVector(7);
theta_7_std = sqrt(model.Report.Parameters.FreeParCovariance(7,7));
disp(['theta_7 = ', num2str(theta_7_mean), ' +/- ', num2str(theta_7_std), ' (V*K)/(A*J)'])

theta_8_mean = model.Report.Parameters.ParVector(8);
theta_8_std = sqrt(model.Report.Parameters.FreeParCovariance(8,8));
disp(['theta_8 = ', num2str(theta_8_mean), ' +/- ', num2str(theta_8_std), '1/s'])

theta_9_mean = model.Report.Parameters.ParVector(9);
theta_9_std = sqrt(model.Report.Parameters.FreeParCovariance(9,9));
disp(['theta_9 = ', num2str(theta_8_mean), ' +/- ', num2str(theta_9_std), '1/s'])

theta_10_mean = model.Report.Parameters.ParVector(10);
theta_10_std = sqrt(model.Report.Parameters.FreeParCovariance(10,10));
disp(['theta_10 = ', num2str(theta_10_mean), ' +/- ', num2str(theta_10_std), '1/s'])

theta_11_mean = model.Report.Parameters.ParVector(11);
theta_11_std = sqrt(model.Report.Parameters.FreeParCovariance(11,11));
disp(['theta_11 = ', num2str(theta_11_mean), ' +/- ', num2str(theta_11_std), '1/s'])

theta_12_mean = model.Report.Parameters.ParVector(12);
theta_12_std = sqrt(model.Report.Parameters.FreeParCovariance(12,12));
disp(['theta_12 = ', num2str(theta_12_mean), ' +/- ', num2str(theta_12_std), '1/s'])

theta_13_mean = model.Report.Parameters.ParVector(13);
theta_13_std = sqrt(model.Report.Parameters.FreeParCovariance(13,13));
disp(['theta_13 = ', num2str(theta_13_mean), ' +/- ', num2str(theta_13_std), '1/s'])

% theta_14_mean = model.Report.Parameters.ParVector(14);
% theta_14_std = sqrt(model.Report.Parameters.FreeParCovariance(14,14));
% disp(['theta_14 = ', num2str(theta_14_mean), ' +/- ', num2str(theta_14_std), '1/s'])

theta_std = [theta_1_std, theta_2_std, theta_3_std, theta_4_std, theta_5_std, ...
    theta_6_std, theta_7_std, theta_8_std, theta_9_std, theta_10_std, theta_11_std, ...
    theta_12_std, theta_13_std];

theta_mean = [theta_1_mean, theta_2_mean, theta_3_mean, theta_4_mean, theta_5_mean, ...
    theta_6_mean, theta_7_mean, theta_8_mean, theta_9_mean, theta_10_mean, ...
    theta_11_mean, theta_12_mean, theta_13_mean];



theta_1 = model.Report.Parameters.ParVector(1);
theta_2 = model.Report.Parameters.ParVector(2);
theta_3 = model.Report.Parameters.ParVector(3);
theta_4 = model.Report.Parameters.ParVector(4);
theta_5 = model.Report.Parameters.ParVector(5);
theta_6 = model.Report.Parameters.ParVector(6);
theta_7 = model.Report.Parameters.ParVector(7);
theta_8 = model.Report.Parameters.ParVector(8);
theta_9 = model.Report.Parameters.ParVector(9);
theta_10 = model.Report.Parameters.ParVector(10);
theta_11 = model.Report.Parameters.ParVector(11);
theta_12 = model.Report.Parameters.ParVector(12);
theta_13 = model.Report.Parameters.ParVector(13);
% theta_14 = model.Report.Parameters.ParVector(14);


%% Abspeichern der Parameter
% ************************************************************************

theta_vektor = [theta_1; theta_2; theta_3; theta_4; theta_5; theta_6; ...
    theta_7; theta_8; theta_9; theta_10; theta_11; theta_12; theta_13];
% 
% save('theta_a3_weighting_open.mat','theta_vektor');



%% Analyse der Parameter (Mean + Standardabweichung):

figure
errorbar((1:13),theta_mean,theta_std)
xlabel('Parameters')
ylabel('Errorbar')
xlim([1,13])


% Mean_Val = mean(theta_vektor, 'all');
% Standard_Deviation_vector = (theta_vektor-Mean_Val).^2;
% 
% figure
% plot((1:numel(theta_vektor)),Standard_Deviation_vector, '.-')
% xlim([1,numel(theta_vektor)])
% ylim([0,max(Standard_Deviation_vector)+max(Standard_Deviation_vector)*0.05])
% xlabel('Parameter')
% ylabel('Standard Deviation')
% title('Standard deviation of each Parameter')