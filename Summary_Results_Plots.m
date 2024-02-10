% This Document delivers a Summary of all results and gives back all
% necessary plots
clear all;
clc;

%% Download of all the relevant files

addpath('Parameter_BSC');
addpath('Parameter_BSC/Messungen_18');
addpath('Theta');

addpath('Data_total');
load('Identification_Response_1.mat');
load('Identification_Response_2.mat');
load('Validation_Response.mat');
load('model_a4_open.mat');
load('Simulation_Result.mat');
load('f_1_a4.mat');

load('c.mat');
alpha_4 = 0.5;
f_alpha_4 = c(2)/alpha_4;

addpath('Data_Door_Opening');
load('Data_CD1_T10.mat');
load('Data_CD1_T14.mat');
load('Data_CD1_T18.mat');
load('Data_CD85_T10.mat');
load('Data_CD85_T14.mat');
load('Data_CD85_T18.mat');
load('par_CD1_T10.mat');
load('par_CD1_T14.mat');
load('par_CD1_T18.mat');
load('par_CD85_T10.mat');
load('par_CD85_T14.mat');
load('par_CD85_T18.mat');

Data_18 = importdata('Experiment_18.mat');
Data_18_new = importdata('Experiment_18_new.mat');

% Note for Identification_1, 2 & Validation:
% Zeile 1: Time
% Zeile 2: Model Response T_HS
% Zeile 3: Measurement Data T_HS
% Zeile 4: Model Response T_ICB
% Zeile 5: Measurement Data T_ICB

% Aufbau sim:
% 1. Spalte: T_wi
% 2. Spalte: T_wo
% 3. Spalte: T_cargo = constant!
% 4. Spalte: T_HS
% 5. Spalte: T_ICB

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
% Aufbau Data_18_new:
% 10. Spalte: P_Fan

% Load of Data not really necessary since I have the Identification and
% Validation Data!
data_1 = load('Messungen_1.mat');
data_1 = data_1.data;
data_2 = load('Messungen_2.mat');
data_2 = data_2.data;
data_3 = load('Messungen_3.mat');
data_3 = data_3.data;

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

% Produce new Data to have every door opening and closing!
door_open_1 = [0];
door_open_2 = [];
door_open_3 = [];

door_close_1 = [];
door_close_2 = [];
door_close_3 = [];

for i = 1:numel(data_1(1,:))
    if i == 1
        prev_val_1 = data_1(3,i);
    else
        if prev_val_1 == data_1(3,i)
            prev_val_1 = data_1(3,i);
        elseif prev_val_1 ~= data_1(3,i)
            if prev_val_1 == 0
                door_open_1 = [door_open_1, data_1(1,i)];
                prev_val_1 = data_1(3,i);
            else
                door_close_1 = [door_close_1, data_1(1,i)];
                prev_val_1 = data_1(3,i);
            end
        end
    end
end

for i = 1:numel(data_2(1,:))
    if i == 1
        prev_val_2 = data_2(3,i);
    else
        if prev_val_2 == data_2(3,i)
            prev_val_2 = data_2(3,i);
        elseif prev_val_2 ~= data_2(3,i)
            if prev_val_2 == 0
                door_open_2 = [door_open_2, data_1(1,i)];
                prev_val_2 = data_2(3,i);
            else
                door_close_2 = [door_close_2, data_2(1,i)];
                prev_val_2 = data_2(3,i);
            end
        end
    end
end

for i = 1:numel(data_3(1,:))
    if i == 1
        prev_val_3 = data_3(3,i);
    else
        if prev_val_3 == data_3(3,i)
            prev_val_3 = data_3(3,i);
        elseif prev_val_3 ~= data_3(3,i)
            if prev_val_3 == 0
                door_open_3 = [door_open_3, data_3(1,i)];
                prev_val_3 = data_3(3,i);
            else
                door_close_3 = [door_close_3, data_3(1,i)];
                prev_val_3 = data_3(3,i);
            end
        end
    end
end

door_open_1 = door_open_1./3600;
door_close_1 = door_close_1./3600;

door_open_2 = door_open_2./3600;
door_close_2 = door_close_2./3600;

door_open_3 = door_open_3./3600;
door_close_3 = door_close_3./3600;

door_open_3_actual = door_open_3 + Identification_1(1,end)/3600 + ...
    Identification_2(1,end)/3600;
door_close_3_actual = door_close_3 + Identification_1(1,end)/3600 + ...
    Identification_2(1,end)/3600;

door_open_tot = [door_open_1, door_open_3_actual];
door_close_tot = [door_close_1, door_close_3_actual];

addpath('Data_Controller');
load('configuration_1.mat');
load('configuration_2.mat');
load('configuration_3.mat');
load('Current.mat');
load('HS_Conf.mat');
load('P_vals_bar.mat');
load('Power_Comparison.mat');
load('Q_z_4_smooth.mat');
load('Residual_Comparison.mat');
load('T_vals_bar.mat');
load('Temp_Violation.mat');
load('door_controller.mat');

% Structure of nearly all data sets just loaded:
% Zeile 1: t
% Zeile 2: Configuration 1
% Zeile 3: Configuration 2
% Zeile 4: Configuration 3


%% Plot Control:
% Turn the values to true if you want to see the plot or store the result!

% Individual Identification Results:
% HW Ratio = 1.3
show_identification_results = false;
store_identification_results = false; % Can only be true if both are true!

% Total Identification Results:
% HW Ratio = 1.3
show_identification_tot_results = false;
store_identification_tot_results = false; % Can only be true if both are true!

% Door Opening:
% NOTE: HWRATIO needs to be adjusted in the function at the end first!
% HW Ratio = 0.85
show_door_opening = false;
store_door_opening = false;

% Regression:
% Note: HW Ratio = 0.5
show_regression = false;
store_regression = false;

% Linearization: 
% HW Ratio = 0.7
show_lin = false;
store_lin = false;

% Parameter Analysis (STD and Evolution!)
% Picture Width: 21/2
show_par_analysis = false;
store_chart_A = false; % HW Ratio = 1.695  (1.72 with 'compact')
store_chart_B = false; % HW Ratio = 1.8

% Controller Configuration Results
show_controller_config = false;
% Residual and HS Plot
% HW Ratio = 0.85
store_Controller_1 = false;
% Plot of all Currents!
% HW Ratio = 0.65
store_Controller_2 = false;

% Controller Configuration 3:
show_config_3 = false;
save_config_3 = false; % HW Ratio: 0.85

% Efficiency Diagram:
% HW Ratio: 1.3
show_efficiency = false;
store_efficiency = false;

% Q_z Diagramm:
% HW Ratio: 1.3
show_Qz = false;
store_Qz = false;

%% Analysis of the Identification Result TOTAL DATASET:

% In this part I want to calculate the Result of just the Identfication
% Datasets and Validation Datasets Seperately!

% This gives back the Analysis Result of the whole Dataset
% 1) ANALYSIS OF HS
% Calculate the R-squared value for each column of data
data_measured_HS = Data_18(:,2);
data_measured_HS = data_measured_HS';
data_estimated_HS = sim(:,4);
data_estimated_HS = data_estimated_HS';


SSres_HS = sum((data_estimated_HS - data_measured_HS).^2);
SStot_HS = sum((data_measured_HS - mean(data_measured_HS)).^2);
R2_HS = 1 - (SSres_HS ./ SStot_HS);

% Display the individual and overall R-squared values
disp('TOTAL DATASET:')
disp(['HS: Individual R-squared values: ', num2str(R2_HS)])

% Calculate NRMSE
range_HS = max(data_measured_HS) - min(data_measured_HS);
nrmse_HS = sqrt(mean((data_measured_HS - data_estimated_HS).^2)) / range_HS;
fit_nrmse_HS = 1-nrmse_HS;
disp(['HS: Overall NRMSE value: ', num2str(fit_nrmse_HS)])

nrmse_HS_mat = 1- sqrt(SSres_HS/SStot_HS);
disp(['HS: Overall NRMSE (Matlab) value: ', num2str(nrmse_HS_mat)])

% 2) ANALYSIS OF ICB
% Calculate the R-squared value for each column of data
data_measured_ICB = Data_18(:,3);
data_measured_ICB = data_measured_ICB';
data_estimated_ICB = sim(:,5);
data_estimated_ICB = data_estimated_ICB';

SSres_ICB = sum((data_measured_ICB - data_estimated_ICB).^2);
SStot_ICB = sum((data_measured_ICB - mean(data_measured_ICB)).^2);
R2_ICB = 1 - (SSres_ICB ./ SStot_ICB);

% Display the individual and overall R-squared values
disp(['ICB: Individual R-squared values: ', num2str(R2_ICB)])

% Calculate NRMSE
range_ICB = max(data_measured_ICB) - min(data_measured_ICB);
nrmse_ICB = sqrt(mean((data_estimated_ICB - data_measured_ICB).^2)) / range_ICB;
fit_nrmse_ICB = 1-nrmse_ICB;
disp(['ICB: Overall NRMSE value: ', num2str(fit_nrmse_ICB)])

nrmse_ICB_mat = 1- sqrt(SSres_ICB/SStot_ICB);
disp(['ICB: Overall NRMSE (Matlab) value: ', num2str(nrmse_ICB_mat)])
disp('_________________________________________________________')

% Clear all irrelevant Data!
clear range_HS nrmse_HS fit_nrmse_HS ...
    data_measured_HS data_estimated_HS SSres_HS SStot_HS ...
    range_ICB nrmse_ICB fit_nrmse_ICB ...
    data_measured_ICB data_estimated_ICB SSres_ICB SStot_ICB

% I need an extra dataset with the splitted results for another plot!
sim_1 = sim(1:numel(Identification_1(1,:)),:);
sim_2 = sim((numel(Identification_1(1,:))+1):2*numel(Identification_1(1,:)),:);
sim_3 = sim((2*numel(Identification_1(1,:))+1):3*numel(Identification_1(1,:)),:);


%% Analysis of the Identification Result Identification Datasets:

% In this part I want to calculate the Result of just the Identfication
% Datasets and Validation Datasets Seperately!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATASET 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) ANALYSIS OF HS
% Calculate the R-squared value for each column of data
data_measured_HS_1 = Identification_1(3,:);
data_estimated_HS_1 = Identification_1(2,:);

SSres_HS_1 = sum((data_estimated_HS_1 - data_measured_HS_1).^2);
SStot_HS_1 = sum((data_measured_HS_1 - mean(data_measured_HS_1)).^2);
R2_HS_1 = 1 - (SSres_HS_1 ./ SStot_HS_1);

% Display the individual and overall R-squared values
disp('IDENTIFICATION DATASET 1:')
disp(['HS: Overall R-squared values: ', num2str(R2_HS_1)])


% Calculate NRMSE
range_HS_1 = max(data_measured_HS_1) - min(data_measured_HS_1);
nrmse_HS_1 = sqrt(mean((data_measured_HS_1 - data_estimated_HS_1).^2)) / ... 
    range_HS_1;
fit_nrmse_HS_1 = 1-nrmse_HS_1;
disp(['HS: Overall NRMSE value: ', num2str(fit_nrmse_HS_1)])

nrmse_HS_mat_1 = 1- sqrt(SSres_HS_1/SStot_HS_1);
disp(['HS: Overall NRMSE (Matlab) value: ', num2str(nrmse_HS_mat_1)])

% 2) ANALYSIS OF ICB
% Calculate the R-squared value for each column of data
data_measured_ICB_1 = Identification_1(5,:);
data_estimated_ICB_1 = Identification_1(4,:);

SSres_ICB_1 = sum((data_estimated_ICB_1 - data_measured_ICB_1).^2);
SStot_ICB_1 = sum((data_measured_ICB_1 - mean(data_measured_ICB_1)).^2);
R2_ICB_1 = 1 - (SSres_ICB_1 ./ SStot_ICB_1);


% Display the individual and overall R-squared values
disp(['ICB: Overall R-squared values: ', num2str(R2_ICB_1)])

% Calculate NRMSE
range_ICB_1 = max(data_measured_ICB_1) - min(data_measured_ICB_1);
nrmse_ICB_1 = sqrt(mean((data_measured_ICB_1 - data_estimated_ICB_1).^2)) /...
    range_ICB_1;
fit_nrmse_ICB_1 = 1-nrmse_ICB_1;
disp(['ICB: Overall NRMSE value: ', num2str(fit_nrmse_ICB_1)])

nrmse_ICB_mat_1 = 1- sqrt(SSres_ICB_1/SStot_ICB_1);
disp(['ICB: Overall NRMSE (Matlab) value: ', num2str(nrmse_ICB_mat_1)])
disp('_________________________________________________________')

% Clear all irrelevant Data!
clear range_HS_1 nrmse_HS_1 fit_nrmse_HS_1 ...
    data_measured_HS_1 data_estimated_HS_1 SSres_HS_1 SStot_HS_1 ...
    range_ICB_1 nrmse_ICB_1 fit_nrmse_ICB_1 ...
    data_measured_ICB_1 data_estimated_ICB_1 SSres_ICB_1 SStot_ICB_1

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATASET 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) ANALYSIS OF HS
% Calculate the R-squared value for each column of data
data_measured_HS_2 = Validation(3,:);
data_estimated_HS_2 = Validation(2,:);

SSres_HS_2 = sum((data_estimated_HS_2 - data_measured_HS_2).^2);
SStot_HS_2 = sum((data_measured_HS_2 - mean(data_measured_HS_2)).^2);
R2_HS_2 = 1 - (SSres_HS_2 ./ SStot_HS_2);

% Display the individual and overall R-squared values
disp('IDENTIFICATION DATASET 2:')
disp(['HS: Overall R-squared values: ', num2str(R2_HS_2)])

% Calculate NRMSE
range_HS_2 = max(data_measured_HS_2) - min(data_measured_HS_2);
nrmse_HS_2 = sqrt(mean((data_measured_HS_2 - data_estimated_HS_2).^2)) / ... 
    range_HS_2;
fit_nrmse_HS_2 = 1-nrmse_HS_2;
disp(['HS: Overall NRMSE value: ', num2str(fit_nrmse_HS_2)])

nrmse_HS_mat_2 = 1- sqrt(SSres_HS_2/SStot_HS_2);
disp(['HS: Overall NRMSE (Matlab) value: ', num2str(nrmse_HS_mat_2)])

% 2) ANALYSIS OF ICB
% Calculate the R-squared value for each column of data
data_measured_ICB_2 = Validation(5,:);
data_estimated_ICB_2 = Validation(4,:);

SSres_ICB_2 = sum((data_estimated_ICB_2 - data_measured_ICB_2).^2);
SStot_ICB_2 = sum((data_measured_ICB_2 - mean(data_measured_ICB_2)).^2);
R2_ICB_2 = 1 - (SSres_ICB_2 ./ SStot_ICB_2);


% Display the individual and overall R-squared values
disp(['ICB: Overall R-squared values: ', num2str(R2_ICB_2)])

% Calculate NRMSE
range_ICB_2 = max(data_measured_ICB_2) - min(data_measured_ICB_2);
nrmse_ICB_2 = sqrt(mean((data_measured_ICB_2 - data_estimated_ICB_2).^2)) /...
    range_ICB_2;
fit_nrmse_ICB_2 = 1-nrmse_ICB_2;
disp(['ICB: Overall NRMSE value: ', num2str(fit_nrmse_ICB_2)])

nrmse_ICB_mat_2 = 1- sqrt(SSres_ICB_2/SStot_ICB_2);
disp(['ICB: Overall NRMSE (Matlab) value: ', num2str(nrmse_ICB_mat_2)])
disp('_________________________________________________________')

% Clear all irrelevant Data!
clear range_HS_2 nrmse_HS_2 fit_nrmse_HS_2 ...
    data_measured_HS_2 data_estimated_HS_2 SSres_HS_2 SStot_HS_2 ...
    range_ICB_2 nrmse_ICB_2 fit_nrmse_ICB_2 ...
    data_measured_ICB_2 data_estimated_ICB_2 SSres_ICB_2 SStot_ICB_2

%%%%%%%%%%%%%%%%%%%%%%%%%% DATASET 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) ANALYSIS OF HS
% Calculate the R-squared value for each column of data
data_measured_HS_3 = Identification_2(3,:);
data_estimated_HS_3 = Identification_2(2,:);

SSres_HS_3 = sum((data_estimated_HS_3 - data_measured_HS_3).^2);
SStot_HS_3 = sum((data_measured_HS_3 - mean(data_measured_HS_3)).^2);
R2_HS_3 = 1 - (SSres_HS_3 ./ SStot_HS_3);

% Display the individual and overall R-squared values
disp('IDENTIFICATION DATASET 3:')
disp(['HS: Overall R-squared values: ', num2str(R2_HS_3)])


% Calculate NRMSE
range_HS_3 = max(data_measured_HS_3) - min(data_measured_HS_3);
nrmse_HS_3 = sqrt(mean((data_measured_HS_3 - data_estimated_HS_3).^2)) / ... 
    range_HS_3;
fit_nrmse_HS_3 = 1-nrmse_HS_3;
disp(['HS: Overall NRMSE value: ', num2str(fit_nrmse_HS_3)])

nrmse_HS_mat_3 = 1- sqrt(SSres_HS_3/SStot_HS_3);
disp(['HS: Overall NRMSE (Matlab) value: ', num2str(nrmse_HS_mat_3)])

% 2) ANALYSIS OF ICB
% Calculate the R-squared value for each column of data
data_measured_ICB_3 = Identification_2(5,:);
data_estimated_ICB_3 = Identification_2(4,:);

SSres_ICB_3 = sum((data_estimated_ICB_3 - data_measured_ICB_3).^2);
SStot_ICB_3 = sum((data_measured_ICB_3 - mean(data_measured_ICB_3)).^2);
R2_ICB_3 = 1 - (SSres_ICB_3 ./ SStot_ICB_3);


% Display the individual and overall R-squared values
disp(['ICB: Overall R-squared values: ', num2str(R2_ICB_3)])

% Calculate NRMSE
range_ICB_3 = max(data_measured_ICB_3) - min(data_measured_ICB_3);
nrmse_ICB_3 = sqrt(mean((data_estimated_ICB_3 - data_measured_ICB_3).^2)) /...
    range_ICB_3;
fit_nrmse_ICB_3 = 1-nrmse_ICB_3;
disp(['ICB: Overall NRMSE value: ', num2str(fit_nrmse_ICB_3)])

nrmse_ICB_mat_3 = 1- sqrt(SSres_ICB_3/SStot_ICB_3);
disp(['ICB: Overall NRMSE (Matlab) value: ', num2str(nrmse_ICB_mat_3)])
disp('_________________________________________________________')

% Clear all irrelevant Data!
clear range_HS_3 nrmse_HS_3 fit_nrmse_HS_3 ...
    data_measured_HS_3 data_estimated_HS_3 SSres_HS_3 SStot_HS_3 ...
    range_ICB_3 nrmse_ICB_3 fit_nrmse_ICB_3 ...
    data_measured_ICB_3 data_estimated_ICB_3 SSres_ICB_3 SStot_ICB_3


% Necesary data for the plots are nrmse_ICB_mat_i 
nrmse_HS_mat = 100 * round(nrmse_HS_mat, 4);
nrmse_ICB_mat = 100 * round(nrmse_ICB_mat, 4);
nrmse_ICB_mat_1 = 100 * round(nrmse_ICB_mat_1, 4);
nrmse_HS_mat_1 = 100 * round(nrmse_HS_mat_1, 4);
nrmse_ICB_mat_2 = 100 * round(nrmse_ICB_mat_2, 4);
nrmse_HS_mat_2 = 100 * round(nrmse_HS_mat_2, 4);
nrmse_ICB_mat_3 = 100 * round(nrmse_ICB_mat_3, 4);
nrmse_HS_mat_3 = 100 * round(nrmse_HS_mat_3, 4);

% Extract Informaiton as String:


text_str_HS = sprintf('Fit: %.2f\\%%', nrmse_HS_mat);
text_str_ICB = sprintf('Fit: %.2f\\%%', nrmse_ICB_mat);

text_str_HS_1 = sprintf('Fit: %.2f\\%%', nrmse_HS_mat_1);
text_str_ICB_1 = sprintf('Fit: %.2f\\%%', nrmse_ICB_mat_1);

text_str_HS_2 = sprintf('Fit: %.2f\\%%', nrmse_HS_mat_2);
text_str_ICB_2 = sprintf('Fit: %.2f\\%%', nrmse_ICB_mat_2);

text_str_HS_3 = sprintf('Fit: %.2f\\%%', nrmse_HS_mat_3);
text_str_ICB_3 = sprintf('Fit: %.2f\\%%', nrmse_ICB_mat_3);

% text_str_HS_3 = "Fit: 78.91\%";
% text_str_ICB_3 = "Fit: 56.06\%";



%% Identification and Validation Dataplot!

% NOTE!!!  set hw_ratio = 1.3

door_color = [255 204 153]./255;
door_alpha = 0.5;

if show_identification_results == true
    
    hfig = figure;
    
    tiledlayout(6,1 ,'TileSpacing','tight','Padding','tight');
    
    nexttile;
    % First Plot Identification 1 and T_HS
    
    if numel(door_open_1) > 0
        for i = 1:numel(door_open_1)
            if i == 1
                x = [door_open_1(i), door_close_1(i), door_close_1(i), ...
                    door_open_1(i)];
                y = [-10, -10, 12, 12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_open_1(i), door_close_1(i), door_close_1(i), ...
                    door_open_1(i)];
                y = [-10, -10, 12, 12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    % Coloured area
    x = [0, Identification_1(1,end), Identification_1(1,end), 0];
    y = [12, 12, 20, 20];
    c = [153 153 255]./255; % Color
    patch(x, y, c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    % Add Text!
    x_center = mean(x(1:2))/3600;
    y_center = mean(y);
    text(x_center, y_center, '$\mathcal{I}_{1}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    p1 = plot(Identification_1(1,:)./3600,Identification_1(2,:)-273.13, 'Color', ...
        [153 153 255]./255,'LineWidth', 2, 'DisplayName', ...
        'Model Response');
    
    p2 = plot(Identification_1(1,:)./3600,Identification_1(3,:)-273.15, 'Color', ...
        [0 0 0],'LineWidth', 1, 'DisplayName', ...
        'Measured Output');
    
    text(2.02, 8, text_str_HS_1, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0]./255, ...
        'FontSize',14)
    
    % TextLocation(text_str_HS_1,'Location', 'best');
    % TextLocation(text_str_HS_1, [0.60, 0.82, 0.1825, 0.0390], 'ident');
    
    hold off
    xlim([0, Identification_1(1,end)/3600])
    ylim([-6, 20])
    y_1 = ylabel('$\vartheta_{cu}$ in $^\circ C$');
    y_1.Position(1) = -0.1;
    label_pos = get(y_1, 'Position');
    
    yticks([-5, 0, 5, 10, 15, 20]); % Set the positions of the ticks
    yticklabels({'' '0' '' '10' '' '20'}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    nexttile
    % Second Plot Identification 1 and T_ICB
    
    if numel(door_open_1) > 0
        for i = 1:numel(door_open_1)
            if i == 1
                x = [door_open_1(i), door_close_1(i), door_close_1(i), ...
                    door_open_1(i)];
                y = [2, 2, 22, 22];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
            else
                x = [door_open_1(i), door_close_1(i), door_close_1(i), ...
                    door_open_1(i)];
                y = [2, 2, 22, 22];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    p1 = plot(Identification_1(1,:)./3600,Identification_1(4,:)-273.13, 'Color', ...
        [153 153 255]./255,'LineWidth', 2, 'DisplayName', ...
        'Model Response');
    
    p2 = plot(Identification_1(1,:)./3600,Identification_1(5,:)-273.15, 'Color', ...
        [0 0 0],'LineWidth', 1, 'DisplayName', ...
        'Measured Output');
    
    text(2.02, 14, text_str_ICB_1, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0]./255, ...
        'FontSize',14)
    
    % TextLocation(text_str_HS_1,'Location', 'best');
    % TextLocation(text_str_HS_1, [0.60, 0.82, 0.1825, 0.0390], 'ident');
    hold off
    xlim([0, Identification_1(1,end)/3600])
    ylim([4, 20.5])
    y_1 = ylabel('$\vartheta_{cc}$ in $^\circ C$');
    y_1.Position(1) = -0.1;
    label_pos = get(y_1, 'Position');
    
    yticks([5, 10, 15, 20]); % Set the positions of the ticks
    yticklabels({'5' '10' '15' '20'}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    
    nexttile
    % Second Identification and T_HS
    if numel(door_open_3) > 0
        for i = 1:numel(door_open_3)
            if i == 1
                x = [door_open_3(i), door_close_3(i), door_close_3(i), ...
                    door_open_3(i)];
                y = [-15, -15, 5.77, 5.77];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_open_3(i), door_close_3(i), door_close_3(i), ...
                    door_open_3(i)];
                y = [-15, -15, 5.77, 5.77];
                c = door_color;
                p = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            end
        end
    end
    x = [0, Identification_1(1,end), Identification_1(1,end), 0];
    y = [5.77, 5.77, 15, 15];
    c = [153 153 255]./255; % Color
    top_2 = patch(x, y, c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    x_center = mean(x(1:2))/3600;
    y_center = mean(y);
    text(x_center, y_center, '$\mathcal{I}_{2}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    
    p3 = plot(Identification_2(1,:)./3600,Identification_2(2,:)-273.15, 'Color', ...
        [153 153 255]./255,'LineWidth', 2, 'DisplayName', ...
        'Model Response');
    hold on
    p4 = plot(Identification_2(1,:)./3600,Identification_2(3,:)-273.15, 'Color', ...
        [0 0 0],'LineWidth', 1, 'DisplayName', ...
        'Measured Output');
    
    text(2.02, 2, text_str_HS_2, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0]./255, ...
        'FontSize',14)
    
    % TextLocation(text_str_HS_3,'Location', 'southwest');
    % TextLocation(text_str_HS_2,[0.1490, 0.45, 0.1825, 0.0390],'ident'); % [0.1490, 0.4227, 0.1825, 0.0390]
    hold off
    xlim([0, Identification_2(1,end)/3600])
    ylim([-15, 15])
    y_2 = ylabel('$\vartheta_{cu}$ in $^\circ C$');
    % y_2.Position(1) = -0.1;
    label_pos_2 = get(y_2, 'Position');
    
    yticks([-15, -10, -5, 0, 5, 10, 15]); % Set the positions of the ticks
    yticklabels({'' '-10' '' '0' '' '10' ''}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    nexttile
    % Second Identification and T_ICB
    if numel(door_open_3) > 0
        for i = 1:numel(door_open_3)
            if i == 1
                x = [door_open_3(i), door_close_3(i), door_close_3(i), ...
                    door_open_3(i)];
                y = [1, 1, 14, 14];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_open_3(i), door_close_3(i), door_close_3(i), ...
                    door_open_3(i)];
                y = [1, 1, 14, 14];
                c = door_color;
                p = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            end
        end
    end
    
    p3 = plot(Identification_2(1,:)./3600,Identification_2(4,:)-273.15, 'Color', ...
        [153 153 255]./255,'LineWidth', 2, 'DisplayName', ...
        'Model Response');
    hold on
    p4 = plot(Identification_2(1,:)./3600,Identification_2(5,:)-273.15, 'Color', ...
        [0 0 0],'LineWidth', 1, 'DisplayName', ...
        'Measured Output');
    
    text(2.02, 8, text_str_ICB_2, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0]./255, ...
        'FontSize',14)
    
    % TextLocation(text_str_HS_3,'Location', 'southwest');
    % TextLocation(text_str_HS_2,[0.1490, 0.45, 0.1825, 0.0390],'ident'); % [0.1490, 0.4227, 0.1825, 0.0390]
    hold off
    xlim([0, Identification_2(1,end)/3600])
    ylim([1, 14])
    y_2 = ylabel('$\vartheta_{cc}$ in $^\circ C$');
    y_2.Position(1) = -0.1;
    label_pos_2 = get(y_2, 'Position');
    
    yticks([2, 4, 6, 8, 10, 12, 14]); % Set the positions of the ticks
    yticklabels({'2' '' '6' '' '10' '' '14'}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    
    nexttile
    % Validation and T_HS
    x = [0, Identification_1(1,end), Identification_1(1,end), 0];
    y = [-1.08, -1.08, 2, 2];
    c = [102,205,170]./255; % Color
    top_3 = patch(x, y, c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on
    x_center = mean(x(1:2))/3600;
    y_center = mean(y);
    text(x_center, y_center, '$\mathcal{V}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    p5 = plot(Validation(1,:)./3600,Validation(2,:)-273.15, 'Color', ...
        [102,205,170]./255,'LineWidth', 2, 'DisplayName', ...
        'Model Response');
    p6 = plot(Validation(1,:)./3600,Validation(3,:)-273.15, 'Color', ...
        [0 0 0],'LineWidth', 1, 'DisplayName', ...
        'Measured Output');
    
    text(2.02, -2.6, text_str_HS_3, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0]./255, ...
        'FontSize',14)
    % TextLocation(text_str_HS_2,'Location', 'best');
    % TextLocation(text_str_HS_3,[0.60, 0.15, 0.1825, 0.0390],'val');
    hold off
    ylim([-8, 2])
    xlim([0, Identification_1(1,end)/3600])
    y_3 = ylabel('$\vartheta_{cu}$ in $^\circ C$');
    y_3.Position(1) = -0.1;
    label_pos_3 = get(y_3, 'Position');
    
    yticks([-8, -6, -4, -2, 0, 2]); % Set the positions of the ticks
    yticklabels({'-8' '' '-4' '' '0' ''}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    nexttile
    % Validation and T_ICB
    p5 = plot(Validation(1,:)./3600,Validation(4,:)-273.15, 'Color', ...
        [102,205,170]./255,'LineWidth', 2, 'DisplayName', ...
        'Model Response');
    hold on
    p6 = plot(Validation(1,:)./3600,Validation(5,:)-273.15, 'Color', ...
        [0 0 0],'LineWidth', 1, 'DisplayName', ...
        'Measured Output');
    
    text(2.02, 4, text_str_ICB_3, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0]./255, ...
        'FontSize',14)
    
    
    % TextLocation(text_str_HS_2,'Location', 'best');
    % TextLocation(text_str_HS_3,[0.60, 0.15, 0.1825, 0.0390],'val');
    hold off
    ylim([-0.5, 6])
    xlim([0, Identification_1(1,end)/3600])
    xlabel('Time in h')
    y_3 = ylabel('$\vartheta_{cc}$ in $^\circ C$');
    y_3.Position(1) = -0.1;
    % label_pos_3 = get(y_3, 'Position');
    
    % Legend
    leg = legend([p1, p5, p2, r],{'Model response of $\mathcal{I}$', ...
        'Model response of $\mathcal{V}$', 'Measured outputs', 'Door opening'}, ...
        'NumColumns',2);
    set(leg,'visible','on');
    set(leg, 'box', 'off');
    set(leg, 'Orientation', 'Vertical');
    leg.Layout.Tile = 'north';
    
    set_figure_properties(hfig);
    
    picturewidth = 21;
    hw_ratio = 1.3; % Figure Height (changeable)
    set(gcf, 'PaperSize', [picturewidth, hw_ratio*picturewidth]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 picturewidth hw_ratio*picturewidth]);
    
    if store_identification_results == true
%         print(gcf, 'pdf_figure', '-dpdf', '-painters', '-fillpage');
        saveas(gcf,'Identification_Individual','pdf');
%         print(hfig, 'pdf_figure', '-dpdf', '-vector', '-fillpage');
    end
    
    clear y_3 y_2 y_1
end

%% Identification Result Total:

door_color = [255 204 153]./255;
door_alpha = 0.5;

if show_identification_tot_results == true
    hfig = figure;
    
    tiledlayout(8,1 ,'TileSpacing','tight','Padding','tight');
    
    nexttile([2,1]);
    % First Plot Identification 1 and T_HS
    if numel(door_open_tot) > 0
        for i = 1:numel(door_open_tot)
            if i == 1
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-15, -15, 15, 15];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-15, -15, 15, 15];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    % Identification Block 1
    x = [0, Identification_1(1,end)/3600, Identification_1(1,end)/3600, 0];
    y = [15, 15, 20.5, 20.5];
    c = [153 153 255]./255; % Color
    patch(x, y, c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    x_center = mean(x(1:2));
    y_center = mean(y);
    text(x_center, y_center, '$\mathcal{I}_{1}$', 'HorizontalAlignment', ...
        'center', 'VerticalAlignment', 'middle')
    
    % Validation Block
    x = [Identification_1(1,end)/3600, 2*Identification_1(1,end)/3600, ...
        2*Identification_1(1,end)/3600, Identification_1(1,end)/3600];
    y = [15, 15, 20.5, 20.5];
    c = [102,205,170]./255; % Color
    patch(x, y, c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    x_center_middle = mean(x(1:2));
    y_center = mean(y);
    text(x_center_middle, y_center, '$\mathcal{V}$', 'HorizontalAlignment', ...
        'center', 'VerticalAlignment', 'middle')
    
    % Identification Block 2
    x = [2*Identification_1(1,end)/3600, 3*Identification_1(1,end)/3600, ...
        3*Identification_1(1,end)/3600, 2*Identification_1(1,end)/3600];
    y = [15, 15, 20.5, 20.5];
    c = [153 153 255]./255; % Color
    patch(x, y, c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    x_center = mean(x(1:2));
    y_center = mean(y);
    text(x_center, y_center, '$\mathcal{I}_{2}$', 'HorizontalAlignment', ...
        'center', 'VerticalAlignment', 'middle')
    
    xline(Identification_1(1,end)/3600, '-.')
    xline(2*Identification_1(1,end)/3600, '-.')
    
    p1 = plot(Identification_1(1,:)./3600,sim_1(:,4), 'Color', ...
        [153 153 255]./255,'LineWidth', 2, 'DisplayName', ...
        'Model Response');
    
    p2 = plot(Data_18((numel(sim_1(:,1))+1):2*numel(sim_1(:,1)),1)/3600, ...
        sim_2(:,4), 'Color', [102,205,170]./255,'LineWidth', 2, ...
        'DisplayName', 'Model Response');
    
    p3 = plot(Data_18((2*numel(sim_1(:,1))+1):3*numel(sim_1(:,1)),1)/3600, ...
        sim_3(:,4), 'Color', [153 153 255]./255,'LineWidth', 2, ...
        'DisplayName', 'Model Response');
    
    p4 = plot(Data_18(:,1)./3600,Data_18(:,2), 'Color', ...
        [0 0 0],'LineWidth', 1, 'DisplayName', ...
        'Measured Output');
    
    text(x_center_middle, 8, text_str_HS, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0], ...
        'FontSize',14)
    
    % TextLocation(text_str_HS_1,'Location', 'best');
    % TextLocation(text_str_HS_1, [0.60, 0.82, 0.1825, 0.0390], 'ident');
    
    hold off
    xlim([0, Data_18(end,1)/3600])
    ylim([-15, 20.5])
    y_1 = ylabel('$\vartheta_{cu}$ in $^\circ C$');
    y_1.Position(1) = -0.4;
    label_pos = get(y_1, 'Position');
    
    % yticks([-5, 0, 5, 10, 15, 20]); % Set the positions of the ticks
    % yticklabels({'' '0' '' '10' '' '20'}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    nexttile([2,1]);
    % Second Plot Identification 1 and T_ICB
    
    if numel(door_open_tot) > 0
        for i = 1:numel(door_open_tot)
            if i == 1
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 22, 22];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 22, 22];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    xline(Identification_1(1,end)/3600, '-.')
    xline(2*Identification_1(1,end)/3600, '-.')
    
    p5 = plot(Identification_1(1,:)./3600,sim_1(:,5), 'Color', ...
        [153 153 255]./255,'LineWidth', 2, 'DisplayName', ...
        'Model Response');
    
    p6 = plot(Data_18((numel(sim_1(:,1))+1):2*numel(sim_1(:,1)),1)/3600, ...
        sim_2(:,5), 'Color', [102,205,170]./255,'LineWidth', 2, ...
        'DisplayName', 'Model Response');
    
    p7 = plot(Data_18((2*numel(sim_1(:,1))+1):3*numel(sim_1(:,1)),1)/3600, ...
        sim_3(:,5), 'Color', [153 153 255]./255,'LineWidth', 2, ...
        'DisplayName', 'Model Response');
    
    p8 = plot(Data_18(:,1)./3600,Data_18(:,3), 'Color', ...
        [0 0 0],'LineWidth', 1, 'DisplayName', ...
        'Measured Output');
    
    text(x_center_middle, 14, text_str_ICB, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0], ...
        'FontSize',14)
    
    % TextLocation(text_str_HS_1,'Location', 'best');
    % TextLocation(text_str_HS_1, [0.60, 0.82, 0.1825, 0.0390], 'ident');
    hold off
    xlim([0, Data_18(end,1)/3600])
    ylim([-1, 20.5])
    y_1 = ylabel('$\vartheta_{cc}$ in $^\circ C$');
    y_1.Position(1) = -0.4;
    label_pos_2 = get(y_1, 'Position');
    
    yticks([0, 5, 10, 15, 20]); % Set the positions of the ticks
    yticklabels({'0' '5' '10' '15' '20'}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    nexttile
    % I Electric Current!
    if numel(door_open_tot) > 0
        for i = 1:numel(door_open_tot)
            if i == 1
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 25, 25];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 25, 25];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    p9 = plot(Data_18(:,1)./3600,Data_18(:,6), 'Color', 'b', ...
        'LineWidth', 1.5, 'DisplayName', 'Model Response');
    
    xline(Identification_1(1,end)/3600, '-.')
    xline(2*Identification_1(1,end)/3600, '-.')
    
    % TextLocation(text_str_HS_3,'Location', 'southwest');
    % TextLocation(text_str_HS_2,[0.1490, 0.45, 0.1825, 0.0390],'ident'); % [0.1490, 0.4227, 0.1825, 0.0390]
    hold off
    xlim([0, Data_18(end,1)/3600])
    ylim([0, 7.5])
    y_2 = ylabel('$I_{tec}$ in A');
    y_2.Position(1) = -0.4;
    label_pos_3 = get(y_2, 'Position');
    
    yticks([0, 2, 4, 6]); % Set the positions of the ticks
    yticklabels({'0' '2' '4' '6'}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    nexttile
    % P_Fan!
    if numel(door_open_tot) > 0
        for i = 1:numel(door_open_tot)
            if i == 1
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 25, 25];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 25, 25];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    p10 = plot(Data_18(:,1)./3600,Data_18_new(:,10), 'Color', 'b', ...
        'LineWidth', 1.5, 'DisplayName', 'Model Response');
    
    xline(Identification_1(1,end)/3600, '-.')
    xline(2*Identification_1(1,end)/3600, '-.')
    
    % TextLocation(text_str_HS_3,'Location', 'southwest');
    % TextLocation(text_str_HS_2,[0.1490, 0.45, 0.1825, 0.0390],'ident'); % [0.1490, 0.4227, 0.1825, 0.0390]
    hold off
    xlim([0, Data_18(end,1)/3600])
    ylim([0, 3.8])
    y_2 = ylabel('$P_{f}$ in W');
    y_2.Position(1) = -0.4;
    label_pos_3 = get(y_2, 'Position');
    
    yticks([0, 1, 2, 3]); % Set the positions of the ticks
    yticklabels({'0' '1' '2' '3'}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    nexttile
    % T_Water!
    if numel(door_open_tot) > 0
        for i = 1:numel(door_open_tot)
            if i == 1
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 25, 25];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 25, 25];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    p11 = plot(Data_18(:,1)./3600,Data_18(:,4), 'Color', [0 0 0], ...
        'LineWidth', 1.5, 'DisplayName', 'Model Response');
    
    xline(Identification_1(1,end)/3600, '-.')
    xline(2*Identification_1(1,end)/3600, '-.')
    
    % TextLocation(text_str_HS_3,'Location', 'southwest');
    % TextLocation(text_str_HS_2,[0.1490, 0.45, 0.1825, 0.0390],'ident'); % [0.1490, 0.4227, 0.1825, 0.0390]
    hold off
    xlim([0, Data_18(end,1)/3600])
    ylim([13, 15.5])
    y_2 = ylabel('$\vartheta_{wtr}$ in $^\circ C$');
    y_2.Position(1) = -0.4;
    label_pos_3 = get(y_2, 'Position');
    
    % yticks([-15, -10, -5, 0, 5, 10, 15]); % Set the positions of the ticks
    % yticklabels({'' '-10' '' '0' '' '10' ''}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    nexttile
    % T_amb
    if numel(door_open_tot) > 0
        for i = 1:numel(door_open_tot)
            if i == 1
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 25, 25];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_open_tot(i), door_close_tot(i), door_close_tot(i),...
                    door_open_tot(i)];
                y = [-10, -10, 25, 25];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    p12 = plot(Data_18(:,1)./3600,Data_18(:,5), 'Color', [0 0 0], ...
        'LineWidth', 1.5, 'DisplayName', 'Model Response');
    
    xline(Identification_1(1,end)/3600, '-.')
    xline(2*Identification_1(1,end)/3600, '-.')
    
    % TextLocation(text_str_HS_3,'Location', 'southwest');
    % TextLocation(text_str_HS_2,[0.1490, 0.45, 0.1825, 0.0390],'ident'); % [0.1490, 0.4227, 0.1825, 0.0390]
    hold off
    xlim([0, Data_18(end,1)/3600])
    ylim([20.5, 22])
    y_2 = ylabel('$\vartheta_{amb}$ in $^\circ C$');
    % y_2.Position(1) = -0.1;
    label_pos_4 = get(y_2, 'Position');
    xlabel('Time in h')
    
    % yticks([2, 4, 6, 8, 10, 12, 14]); % Set the positions of the ticks
    % yticklabels({'2' '' '6' '' '10' '' '14'}); % Set the labels for the ticks
    set(gca, 'TickLength', [0.01 0.01]);
    
    
    
    % Legend
    leg = legend([p1, p2, p4, p9, r],{'Model response of $\mathcal{I}$', ...
        'Model response of $\mathcal{V}$', 'Measured outputs', 'System inputs','Door opening'}, ...
        'NumColumns',3);
    set(leg,'visible','on');
    set(leg, 'box', 'off');
    set(leg, 'Orientation', 'Vertical');
    leg.Layout.Tile = 'north';
    
    set_figure_properties(hfig);
    
    picturewidth = 21;
    hw_ratio = 1.3; % Figure Height (changeable)
    set(gcf, 'PaperSize', [picturewidth, hw_ratio*picturewidth]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 picturewidth hw_ratio*picturewidth]);
    
    
    if store_identification_tot_results == true
        % print(gcf, 'pdf_figure', '-dpdf', '-painters', '-fillpage');
        saveas(gcf,'Identification_Total','pdf');
        % print(hfig, 'pdf_figure', '-dpdf', '-vector', '-fillpage');
    end
    
    clear y_3 y_2 y_1
end

%% Door Opening Plots:

% NOTE!!!  set hw_ratio = 0.85
if show_door_opening == true

    hfig_door = figure;
    
    tiledlayout(2, 2,'TileSpacing', 'compact','Padding','tight');
    
    nexttile
    dummy_large = plot(Data_CD1_T14(3,:), Data_CD1_T14(4,:), "^", ...
        'MarkerEdgeColor', [0, 0, 0]./255);
    hold on
    dummy_small = plot(Data_CD85_T14(3,:), Data_CD85_T14(4,:), "o", ...
        'MarkerEdgeColor', [0, 0, 0]./255);
    p1_1 = plot(Data_CD85_T14(3,:), Data_CD85_T14(4,:), "o", ...
        'MarkerEdgeColor', [26, 137, 237]./255);
    p2_1 = plot(Data_CD1_T14(3,:), Data_CD1_T14(4,:), "^", ...
        'MarkerEdgeColor', [26, 137, 237]./255);
    p1 = plot(Data_CD1_T14(3,:), Data_CD1_T14(4,:), 'LineWidth', 1.5, ...
        'Color', [26, 137, 237]./255);
    p2 = plot(Data_CD85_T14(3,:), Data_CD85_T14(4,:), 'LineWidth', 1.5, ...
        'Color', [26, 137, 237]./255);
    x1_1 = yline(par_CD85_T14.x1_stat, '--', '$x_{1.stat}$', ...
        'LabelHorizontalAlignment', 'center', 'LineWidth', 1.5);
    text(10, 0.28, 'a)', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0], ...
        'FontSize',14)
    hold off
    xlabel('Time in s')
    ylim([0, 0.3])
    y1 = ylabel('Volumetric Flux ($\dot{V}/A_{ap}$) in m/s');
    
    % y_2.Position(1) = -0.4;
    label_pos_1 = get(y1, 'Position');
    
    nexttile
    p3 = plot(Data_CD1_T10(1,:), Data_CD1_T10(2,:), '^', 'MarkerEdgeColor', ...
        [0, 128, 128]./255);
    hold on 
    dummy_3 = plot(Data_CD1_T10(1,:), Data_CD1_T10(2,:), 'Color', ...
        [0, 128, 128]./255, 'LineWidth', 1.5);
    p4 = plot(Data_CD1_T14(1,:), Data_CD1_T14(2,:), '^', 'MarkerEdgeColor', ...
        [26, 137, 237]./255);
    dummy_4 = plot(Data_CD1_T14(1,:), Data_CD1_T14(2,:), 'Color', ...
        [26, 137, 237]./255, 'LineWidth', 1.5);
    p5 = plot(Data_CD1_T18(1,:), Data_CD1_T18(2,:), '^', 'MarkerEdgeColor', ...
        [237, 68, 50]./255);
    dummy_5 = plot(Data_CD1_T18(1,:), Data_CD1_T18(2,:), 'Color', ...
        [237, 68, 50]./255, 'LineWidth', 1.5);
    text(5, 746.67, 'b)', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0], ...
        'FontSize',14)
    hold off
    xlim([0, 75])
    xlabel('Time in s')
    y2 = ylabel('$\dot{Q}_{d}$ in W');
    % y_2.Position(1) = -0.4;
    label_pos_2 = get(y2, 'Position');
    
    nexttile
    % Plot velocity large aperture
    p6 = plot(Data_CD1_T10(3,:), Data_CD1_T10(4,:), '^', 'MarkerEdgeColor', ...
        [0, 128, 128]./255);
    hold on
    p6_1 = plot(Data_CD1_T10(3,:), Data_CD1_T10(4,:), 'Color', ...
        [0, 128, 128]./255, 'LineWidth', 1.5);
    p7 = plot(Data_CD1_T14(3,:), Data_CD1_T14(4,:), '^', 'MarkerEdgeColor', ...
        [26, 137, 237]./255);
    p7_1 = plot(Data_CD1_T14(3,:), Data_CD1_T14(4,:), 'Color', ...
        [26, 137, 237]./255, 'LineWidth', 1.5);
    p8 = plot(Data_CD1_T18(3,:), Data_CD1_T18(4,:), '^', 'MarkerEdgeColor', ...
        [237, 68, 50]./255);
    p8_1 = plot(Data_CD1_T18(3,:), Data_CD1_T18(4,:), 'Color', ...
        [237, 68, 50]./255, 'LineWidth', 1.5);
    text(10, 0.7467, 'c)', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0], ...
        'FontSize',14)
    hold off
    xlabel('Time in s')
    y3 = ylabel('Volumetric Flux ($\dot{V}/A_{ap}$) in m/s');
    y3.Position(1) = -20.5;
    label_pos_3 = get(y3, 'Position');
    
    nexttile
    % Plot Small aperture
    p9 = plot(Data_CD85_T10(3,:), Data_CD85_T10(4,:), 'o', 'MarkerEdgeColor', ...
        [0, 128, 128]./255);
    hold on
    p9_1 = plot(Data_CD85_T10(3,:), Data_CD85_T10(4,:), 'Color', ...
        [0, 128, 128]./255, 'LineWidth', 1.5);
    p10 = plot(Data_CD85_T14(3,:), Data_CD85_T14(4,:), 'o', 'MarkerEdgeColor', ...
        [26, 137, 237]./255);
    p10_1 = plot(Data_CD85_T14(3,:), Data_CD85_T14(4,:), 'Color', ...
        [26, 137, 237]./255, 'LineWidth', 1.5);
    p11 = plot(Data_CD85_T18(3,:), Data_CD85_T18(4,:), 'o', 'MarkerEdgeColor', ...
        [237, 68, 50]./255, 'LineWidth', 1.5);
    p11_1 = plot(Data_CD85_T18(3,:), Data_CD85_T18(4,:), 'Color', ...
        [237, 68, 50]./255, 'LineWidth', 1.5);
    text(10, 0.5133, 'd)', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Color',[0 0 0], ...
        'FontSize',14)
    hold off
    xlabel('Time in s')
    ylim([0, 0.55])
    y4 = ylabel('Volumetric Flux ($\dot{V}/A_{ap}$) in m/s');
    y4.Position(1) = -18;
    label_pos_4 = get(y4, 'Position');
    
    leg = legend([dummy_large, dummy_small, dummy_3, dummy_4, dummy_5], ...
        {'$c_D$ = 1', '$c_D$ = 0.85', '$x_{2.stat}$ = 283.15K', ...
        '$x_{2.stat}$ = 287.15K', '$x_{2.stat}$ = 291.15K'}, ...
        'NumColumns',3);
    set(leg,'visible','on');
    set(leg, 'box', 'on');
    set(leg, 'Orientation', 'Vertical');
    leg.Layout.Tile = 'north';
    
    set_figure_properties(hfig_door);

    if store_door_opening == true
        saveas(gcf,'Door_Opening','pdf');
    end
end


%% Linear Regression Models:

% First Load Data
addpath('Data_System');
% Data to visualize the Linearization!
load('Data_Linearization.mat');
load('I_lin_stat.mat'); % Arbeitspunkt der Linearisierung
% This Data is required for the 3D plot! Those are mesh data!
load('I_grid_sub.mat');
load("Q_3D_sub.mat");
load("f1_fit_sub.mat");
load("T_HS_grid_sub.mat");
% Linear Regression:
load('lin_reg.mat');
load('lin_reg_sub.mat');

% Data_Linearization:
% Data_Linearization(1,:) = I_lin; % x-axis: Current I
% Data_Linearization(2,:) = Q_dot; % Q Values for Linearization
% Data_Linearization(3,:) = Q_dot_lin; % the Linearized Q for a working
% Data_Linearization(4,:) = Data_Q_dot(1,:);
% Data_Linearization(5,:) = Data_Q_dot(2,:);
% Data_Linearization(6,:) = Data_Q_dot(3,:);

% Data Regression:
% lin_reg(1,:) = I_lin;
% lin_reg(2,:) = f2_fit_tot{1};
% lin_reg(3,:) = f2_fit_tot{2};
% lin_reg(4,:) = f2_fit_tot{3};
% lin_reg(5,:) = f3_fit;
% Subdomain of the linear Regression:
% lin_reg_sub(1,:) = I_sub;
% lin_reg_sub(2:4,:) = Q_sub;

% Linearisierung:
% HW = 0.7

if show_lin == true
    idx_find_7A = find(Data_Linearization(1,:) == 7);
    % idx_find_3A = find(Data_Linearization(1,:) == 3);
    
    hfig_lin = figure;
    
    tiledlayout(1, 1,'TileSpacing', 'compact','Padding','tight');
    nexttile
    y1_1 = xline(3, '--','$I_{tec}^{min}$','LabelHorizontalAlignment', ...
        'center','LabelVerticalAlignment', 'bottom','LineWidth', 1.5);
    hold on
    y1_2 = xline(7, '--','$I_{tec}^{max}$','LabelHorizontalAlignment', ...
        'center','LabelVerticalAlignment', 'bottom','LineWidth', 1.5);
    
    p2 = xline(I_lin_stat,'-', '$I_{tec}^{lin}$','LabelHorizontalAlignment',...
        'center', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
    p1 = plot(Data_Linearization(1,:) , Data_Linearization(2,:), ...
        'Color',[135 159 255]./255, 'LineWidth', 2);
    p4 = plot(Data_Linearization(1, 1:idx_find_7A), ...
        Data_Linearization(3, 1:idx_find_7A),'Color', [107 153 50]./255, ...
        'LineWidth', 2);
    hold off
    xlim([0,8])
    ylabel('$\dot{Q}$ in kW', 'Interpreter', 'latex')
    xlabel('$I_{tec}$ in A', 'Interpreter', 'latex')
    leg = legend([p1, p4],{'$\dot{Q}_{sub}^{nonlin}$', '$\dot{Q}_{sub}^{lin}$'},...
        'Interpreter','latex');
    set(leg,'visible','on');
    set(leg, 'box', 'on');
    set(leg, 'Orientation', 'Horizontal');
    leg.Layout.Tile = 'north';
    
    set_figure_properties(hfig_lin);


    if store_lin == true
        saveas(gcf,'Linearization','pdf');
    end
end



% Regression Models in 1 image!
% HW = 0.5
if show_regression == true
    idx_find_7A = find(lin_reg(1,:) == 7);
    idx_find_3A = find(lin_reg(1,:) == 3);
    
    
    hfig_reg = figure;
    
    tiledlayout(1, 2,'TileSpacing', 'loose','Padding','tight');
    
    nexttile
    p1 = surf(I_grid_sub, T_HS_grid_sub, Q_3D_sub,'FaceColor',...
        [135 159 255]./255, 'EdgeColor', 'none'); % 230 109 96
    hold on
    p2 = surf(I_grid_sub, T_HS_grid_sub, f1_fit,'FaceColor', ...
        [230 109 96]./255, 'EdgeColor', 'none');
    hold off
    xticks([3 5 7])
    % legend('Original Data', 'Fitted Plane')
    hl = xlabel('$I_{tec}$ in A'); % label x-axis              
    vl = ylabel('$T_{cu.stat}$ in K'); % label y-axis
    set(hl, 'Position', get(hl, 'Position') - [0 0.1 0])
    set(vl, 'Position', get(vl, 'Position') - [0.1 0 0])
    zlabel('$\dot{Q}$ in kW') % label z-axis
    
    nexttile
    p3 = plot(Data_Linearization(1,:), Data_Linearization(4,:),'Color', ...
        [135 159 255]./255, 'LineWidth', 1.0);
    hold on;
    p3_1 = plot(Data_Linearization(1,:), Data_Linearization(5,:),'Color', ...
        [135 159 255]./255, 'LineWidth', 1.0);
    p3_2 = plot(Data_Linearization(1,:), Data_Linearization(6,:),'Color', ...
        [135 159 255]./255, 'LineWidth', 1.0);
    % hold on;
    p5 = plot(lin_reg(1, idx_find_3A:idx_find_7A), ...
        lin_reg(2, idx_find_3A:idx_find_7A),'Color', ...
        [230 109 96]./255, 'LineWidth', 1.5);
    p5_1 = plot(lin_reg(1, 1:idx_find_3A), ...
        lin_reg(2, 1:idx_find_3A), '--','Color', ...
        [230 109 96]./255, 'LineWidth', 1.0);
    p6 = plot(lin_reg(1, idx_find_3A:idx_find_7A), ...
        lin_reg(3, idx_find_3A:idx_find_7A),'Color', ...
        [230 109 96]./255, 'LineWidth', 1.5);
    p6_1 = plot(lin_reg(1, 1:idx_find_3A), ...
        lin_reg(3, 1:idx_find_3A), '--','Color', ...
        [230 109 96]./255, 'LineWidth', 1.0);
    p7 = plot(lin_reg(1, idx_find_3A:idx_find_7A), ...
        lin_reg(4, idx_find_3A:idx_find_7A),'Color', ...
        [230 109 96]./255, 'LineWidth', 1.5);
    p7_1 = plot(lin_reg(1, 1:idx_find_3A), ...
        lin_reg(4, 1:idx_find_3A), '--','Color', ...
        [230 109 96]./255, 'LineWidth', 1.0);
    p4 = plot(lin_reg(1,1:idx_find_7A), lin_reg(5,1:idx_find_7A),'Color', ...
        [107 153 50]./255, 'LineWidth', 1.5);
    y1_1 = xline(3, '--', '$I_{tec}^{min}$', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', ...
        'center', 'LineWidth', 1.5);
    y1_2 = xline(7, '--', '$I_{tec}^{max}$', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', ...
        'center', 'LineWidth', 1.5);
    hold off
    xlim([0, 8])
    ylim([0, 0.155])
    xl = xlabel({'$I_{tec}$ in A';'b) Additional linearization methods'});
    ylabel('$\dot{Q}$ in kW')
    pos_xl = get(xl, 'Position');
    leg = legend([p3, p7, p4], ...
        {'$\dot{Q}_{sub}$', '$\dot{Q}_{sub}^{lin}$(.;1)', ...
        '$\dot{Q}_{sub}^{lin}$(.;0)'});
    set(leg,'visible','on');
    set(leg, 'box', 'on');
    set(leg, 'Orientation', 'Horizontal');
    leg.Layout.Tile = 'north';
    
    % Add text annotation to first plot
    ann_xpos = 0.25;
    ann_ypos = pos_xl(2);
    annotation('textbox', [0.2 0.04 0.1 0.1], ...
        'String', {'';'';'a) Linearization method in space'}, ...
        'HorizontalAlignment', 'center', ...
        'Color', [0 0 0], ...
        'EdgeColor', 'none')
    
    set_figure_properties(hfig_reg);

    if store_regression == true
        saveas(gcf,'Regression_Models','pdf');
    end
end

%% Parameter Estimation:

if show_par_analysis == true

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
    
    theta_std = [theta_1_std, theta_2_std, theta_3_std, theta_4_std, theta_5_std, ...
        theta_6_std, theta_7_std, theta_8_std, theta_9_std, theta_10_std, theta_11_std, ...
        theta_12_std, theta_13_std];
    
    theta_mean = [theta_1_mean, theta_2_mean, theta_3_mean, theta_4_mean, theta_5_mean, ...
        theta_6_mean, theta_7_mean, theta_8_mean, theta_9_mean, theta_10_mean, ...
        theta_11_mean, theta_12_mean, theta_13_mean];
    
    theta_1_std_percent = 100 * theta_1_std / theta_1_mean;
    theta_2_std_percent = 100 * theta_2_std / theta_2_mean;
    theta_3_std_percent = 100 * theta_3_std / theta_3_mean;
    theta_4_std_percent = 100 * theta_4_std / theta_4_mean;
    theta_5_std_percent = 100 * theta_5_std / theta_5_mean;
    theta_6_std_percent = 100 * theta_6_std / theta_6_mean;
    theta_7_std_percent = 100 * theta_7_std / theta_7_mean;
    theta_8_std_percent = 100 * theta_8_std / theta_8_mean;
    theta_9_std_percent = 100 * theta_9_std / theta_9_mean;
    theta_10_std_percent = 100 * theta_10_std / theta_10_mean;
    theta_11_std_percent = 100 * theta_11_std / theta_11_mean;
    theta_12_std_percent = 100 * theta_12_std / theta_12_mean;
    theta_13_std_percent = 100 * theta_13_std / theta_13_mean;
    
    theta_std_percent = [theta_1_std_percent, theta_2_std_percent, ...
        theta_3_std_percent, 0, theta_4_std_percent, theta_5_std_percent, ...
        theta_6_std_percent, theta_7_std_percent, theta_8_std_percent, ...
        theta_9_std_percent, theta_10_std_percent, theta_11_std_percent, ...
        theta_12_std_percent, 0];
    
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
    
    theta_vektor = [theta_1; theta_2; theta_3; theta_4; theta_5; theta_6; ...
        theta_7; theta_8; theta_9; theta_10; theta_11; theta_12; theta_13];
    th = theta_vektor; % Normalized
    
    th_tot = theta_vektor.*f_1;
    th_all = [th_tot(1), th_tot(2), th_tot(3), c(2), th_tot(4), th_tot(5), ...
        th_tot(6), th_tot(7), th_tot(8), th_tot(9), th_tot(10), th_tot(11), ...
        th_tot(12), th_tot(13)]; % Non normalized results!
    
    theta_tot = [th(1), th(2), th(3), alpha_4, th(4), th(5), th(6), th(7), ...
        th(8), th(9), th(10), th(11), th(12), th(13)]; % normalized results
    
    theta_tot_old = ones(1,14).*0.5; % Normalized Values OLD (0.5)
    clear th;
    
    th_old = 0.5.*f_1;
    th_old_all = [th_old(1), th_old(2), th_old(3), c(2), th_old(4), ...
        th_old(5), th_old(6), th_old(7), th_old(8), th_old(9), th_old(10), ...
        th_old(11), th_old(12), th_old(13)]; % Non normalized values OLD!
    
    th_analysis(1,:) = th_all; % Non normalized Values NEW
    th_analysis(2,:) = th_old_all; % Non normalized Values OLD
    
    % Again without alpha_3!!! Due to outlier Problem!
    idx_plot = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
    idx_helper = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
    
    th_analysis_helper = th_analysis;
    th_analysis(:,4) = -10;
    th_analysis_helper(:, idx_helper) = -10;
    th_diff = th_analysis(2,idx_plot)-th_analysis(1,idx_plot);
    th_helper(:,1) = th_analysis(1,idx_plot);
    th_helper(:,2) = th_diff;
    
    % HW Ratio = 1.695
    chart_1 = figure;
    t_1 = tiledlayout(1, 1, 'Padding', 'compact');
    ax1 = axes(t_1);
    yline(1:13, 'LineWidth',0.25,'color',[.75 .75 .75])
    hold on
    h = barh(ax1, 1:13, th_helper, 'stacked');
    h(2).FaceColor = [.75 .75 .75];
    h(2).FaceAlpha =  0.3;
    h(2).BarWidth = 0.3;
    h(2).EdgeColor = 'none';
    h(1).Visible = 'off';
    p2 = plot(ax1, th_analysis(1,idx_plot), 1:13, 'o', 'MarkerSize', 10, ...
        "MarkerEdgeColor","none","MarkerFaceColor", [255 145 77]./255);
    p1 = plot(ax1, th_analysis(2,idx_plot), 1:13, 'o', 'MarkerSize', 10, ...
        "MarkerEdgeColor","none","MarkerFaceColor", [204 78 0]./255);
    p3 = plot(ax1, th_analysis_helper(1, idx_plot), 1:13, 'diamond',...
        'MarkerSize', 10, "MarkerEdgeColor","none","MarkerFaceColor", ...
        [255 0 0]./255);
    hold off
    xlabel('Parameter value in 1')
    ylabel('Parameter')
    title({'\rm a) Comparison of the parameter values'; ''})
    xlim([0, 0.21])
    ylim([0, 13.5])
    yticks(1:13)
    yticklabels(["$\alpha_1$", "$\alpha_2$", "$\alpha_4$", "$\beta_1$", ...
        "$\beta_2$","$\beta_3$","$\beta_4$","$\beta_5$","$\gamma_1$",...
        "$\gamma_2$", "$\gamma_3$", "$\gamma_4$", "$\gamma_5$",])
    leg = legend([p1, p2, p3], ...
        {'Initial value', 'Estimated value', 'Fixed value'});
    set(leg,'visible','on');
    set(leg, 'box', 'off');
    set(leg, 'Orientation', 'Vertical');
    leg.Layout.Tile = 'south';
    
    ax2 = axes(t_1);
    xlim([0, 0.21])
    ylim([0, 13.5])
    yticks(1:13)
    % set(gca,'YMinorTick','on')
    set(gca,'yticklabel',{[]})
    set(gca,'xticklabel',{[]})
    
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
    ax2.Color = 'none';
    ax1.Box = 'off';
    ax2.Box = 'off';
    
    set_figure_properties(chart_1);
    if store_chart_A == true
        saveas(gcf,'Barchart_A','pdf');
    end
    
    % Relevant Idx for the plot:
    std_helper = zeros(1, 13);
    th_tot_helper_1 = theta_tot;
    th_tot_helper_1(4) = -10;
    
    th_tot_helper_2 = theta_tot;
    th_tot_helper_2(idx_plot) = -10;
    th_tot_helper_2(4) = theta_tot(4);
    
    % HW Ratio = 1.8
    chart_2 = figure;
    t_2 = tiledlayout(1, 1, 'Padding', 'compact');
    ax1 = axes(t_2);
    yline(1:13, 'LineWidth',0.25,'color',[.75 .75 .75])
    hold on
    p1 = barh(ax1, 1:13, theta_std_percent(idx_plot).*(0.6/125), ...
        'EdgeColor', 'none', 'FaceColor', "#80B3FF", 'FaceAlpha', 0.3);
    yline(1:13, 'LineWidth',0.25,'color',[.75 .75 .75])
    p2 = plot(ax1, th_tot_helper_1(idx_plot) , 1:13, "o", 'MarkerSize', 10,...
        "MarkerEdgeColor","none","MarkerFaceColor", [79 120 255]./255);
    % [0.65 0.85 0.90]
    hold on
    p3 = plot(ax1, th_tot_helper_2(idx_plot) , 1:13, "diamond",...
        'MarkerSize', 10, "MarkerEdgeColor","none","MarkerFaceColor", [0 0 1]);
    hold off
    xlabel('Estimated normalized value in 1', ...
        'Interpreter', 'latex') % 'b) Estimated values and standard deviation'
    ylabel('Normalized parameter')
    xlim([0, 0.6])
    ylim([0, 13.5])
    xticks([0 0.1 0.2 0.3 0.4 0.5 0.6])
    yticks(1:13)
    yticklabels(["$\alpha_1'$", "$\alpha_2'$", "$\alpha_4'$", "$\beta_1'$", ...
        "$\beta_2'$","$\beta_3'$","$\beta_4'$","$\beta_5'$","$\gamma_1'$",...
        "$\gamma_2'$", "$\gamma_3'$", "$\gamma_4'$", "$\gamma_5'$",])
    
    leg = legend([p1, p2, p3], ...
        {'Standard deviation', 'Estimated value', 'Fixed value'});
    set(leg,'visible','on');
    set(leg, 'box', 'on');
    set(leg, 'Orientation', 'Vertical');
    leg.Layout.Tile = 'south';
    
    ax2 = axes(t_2);
    xlim([0, 125])
    ylim([0, 13.5])
    yticks(1:13)
    % set(gca,'YMinorTick','on')
    xticks([0, 20, 40, 60, 80, 100, 120])
    xlabel('Estimated standard deviation in $\%$ of mean', 'Interpreter', ...
        'latex')
    set(gca,'yticklabel',{[]})
    title('\rm b) Estimated values and standard deviation')
    
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
    ax2.Color = 'none';
    ax1.Box = 'off';
    ax2.Box = 'off';
    set_figure_properties(chart_2);
    
    if store_chart_B == true
        saveas(gcf,'Barchart_B','pdf');
    end
end


%% Controller Description

door_controller = door_controller ./3600;

if show_controller_config == true

    % Two figures required:
    %HW Ratio = 0.85
  
    
        
    hfig_conf = figure;
    
    tiledlayout(2,1 ,'TileSpacing','compact','Padding','tight');
    
    nexttile;
    % First Plot: Residual Plot with all Configurations!
    
    
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
                x_hyst = [0, door_controller(1,i), door_controller(1,i), 0];
                y_hyst = [-0.5, -0.5, 0.5, 0.5];
                c_hyst = [.5 .5 .5]; % Color
                patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
                x_hyst = [door_controller(2,i), door_controller(1,i+1), ...
                    door_controller(1,i+1), door_controller(2,i)];
                hyst = patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none');
    
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                
                if i < numel(door_controller(1,:))
                    x_hyst = [door_controller(2,i), door_controller(1,i+1), ...
                        door_controller(1,i+1), door_controller(2,i)];
                else
                    x_hyst = [door_controller(2,i), 8, 8, door_controller(2,i)];
                end
    
                y_hyst = [-0.5, -0.5, 0.5, 0.5];
                c_hyst = [.5 .5 .5]; % Color
                patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    
            end
        end
    end
    
    x_1 = [0, 2, 2, 0];
    x_2 = [2, 4, 4, 2];
    x_3 = [4, 6, 6, 4];
    x_4 = [6, 8, 8, 6];
    y = [13.12, 13.12, 15, 15];
    
    c_1 = [15, 194, 192]./255; % Color
    c_2 = [0, 143, 140]./255; % Color
    c_3 = [1, 89, 88]./255; % Color
    c_4 = [2, 53, 53]./255; % Color
    
    patch(x_1, y, c_1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_2, y, c_2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_3, y, c_3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_4, y, c_4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Add Text!
    x_center_1 = mean(x_1(1:2));
    x_center_2 = mean(x_2(1:2));
    x_center_3 = mean(x_3(1:2));
    x_center_4 = mean(x_4(1:2));
    y_center = mean(y);
    text(x_center_1, y_center, '$\mathcal{S}_{1}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_2, y_center, '$\mathcal{S}_{2}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_3, y_center, '$\mathcal{S}_{3}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_4, y_center, '$\mathcal{S}_{4}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    
    xline(2, '-.')
    xline(4, '-.')
    xline(6, '-.')
    
    p3 = plot(Residual_Comparison(1,:), Residual_Comparison(4,:), 'Color', ...
        '#d95555','LineWidth', 1.5, 'DisplayName', ...
        'Measured Output');
    
    p2 = plot(Residual_Comparison(1,:), Residual_Comparison(3,:), 'Color', ...
        '#7EB283','LineWidth', 1.5, 'DisplayName', ...
        'Measured Output'); 
    
    p1 = plot(Residual_Comparison(1,:), Residual_Comparison(2,:), 'Color', ...
        '#1a46a3','LineWidth', 1.5, 'DisplayName', ...
        'Model Response');
    
    
    
    hold off
    % set(gca,'YMinorTick','on')
    % set(gca,'XMinorTick','on')
    % set(gca, 'TickLength', [0.01 0.01]);
    xlim([0, 8])
    ylim([-1, 15])
    y_1 = ylabel('$r_{cc}$ in $^\circ C$');
    
    nexttile
    % Second Plot
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 20, 20];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 20, 20];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    xline(2, '-.')
    xline(4, '-.')
    xline(6, '-.')
    
    p6 = plot(HS_Conf(1,:), HS_Conf(4,:), 'Color', ...
        '#d95555','LineWidth', 1.5, 'DisplayName', ...
        'Measured Output');
    
    p5 = plot(HS_Conf(1,:), HS_Conf(3,:), 'Color', ...
        '#7EB283','LineWidth', 1.5, 'DisplayName', ...
        'Measured Output');
    
    p4 = plot(HS_Conf(1,:), HS_Conf(2,:), 'Color', ...
        '#1a46a3','LineWidth', 1.5, 'DisplayName', ...
        'Model Response'); 
    
    hold off
    xlim([0, 8])
    ylim([-8, 20])
    y_2 = ylabel('$\vartheta_{cu}$ in $^\circ C$');
    xlabel('Time in h')
    
    % set(gca, 'TickLength', [0.01 0.01]);
    
    leg = legend([p1, p2, p3, hyst, r], ...
        {'Controller $\mathcal{C}_{1}$', 'Controller $\mathcal{C}_{2}$', ...
        'Controller $\mathcal{C}_{3}$', 'Hysteresis band', 'Door opening'}, ...
        'NumColumns',3);
    set(leg,'visible','on');
    set(leg, 'box', 'off');
    set(leg, 'Orientation', 'Vertical');
    leg.Layout.Tile = 'north';
        
    set_figure_properties(hfig_conf)
    
    if store_Controller_1 == true
        saveas(gcf,'Residual_HS_Controller','pdf');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Second Plot:
    % All Currents of all Configurations!
    % HW Ratio = 0.65
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    hfig_conf_2 = figure;
    
    tiledlayout(3,1 ,'TileSpacing','compact','Padding','tight');
    
    nexttile;
    % Current I of Controller Configuration 1
    
    
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    p1 = plot(Current(1,:), Current(2,:), 'Color', 'blue', 'LineWidth', 1);
    
    xlim([0, 8])
    ylim([0, 7.5])
    yticks([0, 2, 4, 6])
    title('\rm a) Controller $\mathcal{C}_{1}$')
    ylabel('$I_{tec}$ in A')
    
    nexttile
    
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    p2 = plot(Current(1,:), Current(3,:), 'Color', 'blue', 'LineWidth', 1);
    
    xlim([0, 8])
    ylim([0, 7.5])
    yticks([0, 2, 4, 6])
    title('\rm b) Controller $\mathcal{C}_{2}$')
    ylabel('$I_{tec}$ in A')
    
    nexttile
    
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    p3 = plot(Current(1,:), Current(4,:), 'Color', 'blue', 'LineWidth', 1);
    
    xlim([0, 8])
    ylim([0, 7.5])
    yticks([0, 2, 4, 6])
    title('\rm c) Controller $\mathcal{C}_{3}$')
    ylabel('$I_{tec}$ in A')
    xlabel('Time in h')
    
    leg = legend([p1, r], ...
        {'System input', 'Door opening'});
    set(leg,'visible','on');
    set(leg, 'box', 'off');
    set(leg, 'Orientation', 'Horizontal');
    leg.Layout.Tile = 'north';
    
    set_figure_properties(hfig_conf_2)
    
    if store_Controller_2 == true
        saveas(gcf,'Current_I_Outputs_2','pdf');
    end
end

%% Relevant Configuration 3 Plots

% HW Ratio: 0.85
if show_config_3 == true
    idx_t = find(Current(1,:) == 1.1);
    
    TRef = 4.5; %Set reference
    TDiff = abs(configuration_3.T_ICB.Data(1:idx_t) - TRef);
    [minVal, minInd] = min(TDiff); %Find value closest to TRef
    
    
    hfig_conf_3 = figure;
    
    tiledlayout(5,1 ,'TileSpacing','compact','Padding','tight');
    
    nexttile([2 1]);
        
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12+5, 13.12+5];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
                x_hyst = [0, door_controller(1,i), door_controller(1,i), 0];
                y_hyst = [4.5, 4.5, 5.5, 5.5];
                c_hyst = [.5 .5 .5]; % Color
                patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
                x_hyst = [door_controller(2,i), door_controller(1,i+1), ...
                    door_controller(1,i+1), door_controller(2,i)];
                hyst = patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none');
    
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12+5, 13.12+5];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
                if i < numel(door_controller(1,:))
                    x_hyst = [door_controller(2,i), door_controller(1,i+1), ...
                        door_controller(1,i+1), door_controller(2,i)];
                else
                    x_hyst = [door_controller(2,i), 8, 8, door_controller(2,i)];
                end
    
                y_hyst = [4.5, 4.5, 5.5, 5.5];
                c_hyst = [.5 .5 .5]; % Color
                patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    
            end
        end
    end
    
    x_1 = [0, 2, 2, 0];
    x_2 = [2, 4, 4, 2];
    x_3 = [4, 6, 6, 4];
    x_4 = [6, 8, 8, 6];
    y = [13.12+5, 13.12+5, 15+5, 15+5];
    
    c_1 = [15, 194, 192]./255; % Color
    c_2 = [0, 143, 140]./255; % Color
    c_3 = [1, 89, 88]./255; % Color
    c_4 = [2, 53, 53]./255; % Color
    
    patch(x_1, y, c_1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_2, y, c_2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_3, y, c_3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_4, y, c_4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Add Text!
    x_center_1 = mean(x_1(1:2));
    x_center_2 = mean(x_2(1:2));
    x_center_3 = mean(x_3(1:2));
    x_center_4 = mean(x_4(1:2));
    y_center = mean(y);
    
    text(x_center_1, y_center, '$\mathcal{S}_{1}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_2, y_center, '$\mathcal{S}_{2}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_3, y_center, '$\mathcal{S}_{3}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_4, y_center, '$\mathcal{S}_{4}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    
    xline(2, '-.')
    xline(4, '-.')
    xline(6, '-.')
    % line([minInd minInd],[0, 13.12+5],'Color','k');
    
    % x_problem = [Current(1,minInd), 1.8833, 1.8833, Current(1,minInd)];
    % y_problem = [5.5, 5.5, 13.12+5, 13.12+5];
    % c_problem = [2, 53, 53]./255; % Color
    % patch(x_problem, y_problem, c_problem, 'FaceAlpha', 0.2, 'EdgeColor', ...
    %     'none');
    
    p1 = plot(Current(1,:), configuration_3.T_ICB.Data, 'Color', '#A60A33', ...
        'LineWidth', 1.5); % #591E29
    p2 = plot(Current(1,:), configuration_3.T_wi.Data, 'Color', '#D9ADAD', ...
        'LineWidth', 1.5);
    
    dim = [.178 .601 .03 .184];
    annotation('ellipse',dim)
    
    hold off
    
    xlim([0 8])
    ylim([4, 20])
    
    y1 = ylabel('$\vartheta$ in $^\circ C$');
    label_pos_1 = get(y1, 'Position');
    
    
    nexttile([2 1]);
    
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    xline(2, '-.')
    xline(4, '-.')
    xline(6, '-.')
    
    p3 = plot(Current(1,:), configuration_3.PI.Data, 'Color',[0 0 0], ...
        'LineWidth', 1.5);
    p4 = plot(Current(1,:), Current(4,:), 'Color', 'blue', 'LineWidth', 1.5);
    
    hold off
    
    ylim([0, 7.5])
    xlim([0 8])
    
    y2 = ylabel('$I_{tec}$ in A');
    label_pos_2 = get(y2, 'Position');
    y2.Position(1) = -0.3;
    
    nexttile
    
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [-10, -10, 13.12, 13.12];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    
    xline(2, '-.')
    xline(4, '-.')
    xline(6, '-.')
    
    p5 = plot(Current(1,:), configuration_3.s_fan.Data, 'Color', 'blue', ...
        'LineWidth', 1.5);
    
    hold off
    yticks([0 1])
    ylim([0, 1.1])
    xlim([0 8])
    
    xlabel('Time in h')
    y3 = ylabel('$s_{f}$ in 1');
    label_pos_3 = get(y3, 'Position');
    y3.Position(1) = -0.3;
    
    leg = legend([p1, p2, p3, p4, hyst, r], ...
        {'$\vartheta_{cc}$', '$\vartheta_{w.1}$','PI controller output', ...
        'System inputs', 'Hysteresis band', 'Door opening'}, 'NumColumns', 3);
    set(leg,'visible','on');
    set(leg, 'box', 'off');
    set(leg, 'Orientation', 'Vertical');
    leg.Layout.Tile = 'north';
    
    set_figure_properties(hfig_conf_3)
    
    if save_config_3 == true
        saveas(gcf,'Configuration_3','pdf');
    end
end




%% Energy Efficiency

if store_efficiency == true
    
    hfig_efficiency = figure;
    
    tiledlayout(3,2 ,'TileSpacing','compact','Padding','tight');
    
    nexttile
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [0, -0,  4.23077, 4.23077].*(10^5);
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [0, -0,  4.23077, 4.23077].*(10^5);
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    x_1 = [0, 2, 2, 0];
    x_2 = [2, 4, 4, 2];
    x_3 = [4, 6, 6, 4];
    x_4 = [6, 8, 8, 6];
    y = [4.23077, 4.23077, 5, 5].*(10^5);
    
    c_1 = [15, 194, 192]./255; % Color
    c_2 = [0, 143, 140]./255; % Color
    c_3 = [1, 89, 88]./255; % Color
    c_4 = [2, 53, 53]./255; % Color
    
    patch(x_1, y, c_1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_2, y, c_2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_3, y, c_3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_4, y, c_4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Add Text!
    x_center_1 = mean(x_1(1:2));
    x_center_2 = mean(x_2(1:2));
    x_center_3 = mean(x_3(1:2));
    x_center_4 = mean(x_4(1:2));
    y_center = mean(y);
    
    text(x_center_1, y_center, '$\mathcal{S}_{1}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_2, y_center, '$\mathcal{S}_{2}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_3, y_center, '$\mathcal{S}_{3}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_4, y_center, '$\mathcal{S}_{4}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    
    xline(2, '-.')
    xline(4, '-.')
    xline(6, '-.')
    
    % Power_Comparison(1,:) = t;
    % Power_Comparison(2,:) = P_I_1;
    % Power_Comparison(3,:) = P_I_2;
    % Power_Comparison(4,:) = P_I_3;
    
    p1 = plot(Power_Comparison(1,:), Power_Comparison(2,:), 'Color', ...
        '#d95555', 'LineWidth', 1.5);
    
    p2 = plot(Power_Comparison(1,:), Power_Comparison(3,:), 'Color', ...
        '#7EB283', 'LineWidth', 1.5);
    
    p3 = plot(Power_Comparison(1,:), Power_Comparison(4,:), 'Color', ...
        '#1a46a3', 'LineWidth', 1.5);
    hold off
    
    P_max_1 = max(Power_Comparison(2,:));
    P_max_2 = max(Power_Comparison(3,:));
    P_max_3 = max(Power_Comparison(4,:));
    P_max = max([P_max_1, P_max_2, P_max_3]);
    
    xlim([0, 8])
    ylim([0, 5*10^5])
    
    xlabel('Time in h')
    y1 = ylabel('$\hat{P}$ in $A^2$s');
    
    y1.Position(1) = -0.8;
    label_pos_1 = get(y1, 'Position');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    if numel(door_open_1) > 0
        for i = 1:numel(door_controller(1,:))
            if i == 1
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [0, 0,  6.7692, 6.7692];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                hold on
    
            else
                x = [door_controller(1,i), door_controller(2,i), door_controller(2,i), ...
                    door_controller(1,i)];
                y = [0, 0,  6.7692, 6.7692];
                c = door_color;
                r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
            end
        end
    end
    
    x_1 = [0, 2, 2, 0];
    x_2 = [2, 4, 4, 2];
    x_3 = [4, 6, 6, 4];
    x_4 = [6, 8, 8, 6];
    y = [6.7692, 6.7692, 8, 8];
    
    c_1 = [15, 194, 192]./255; % Color
    c_2 = [0, 143, 140]./255; % Color
    c_3 = [1, 89, 88]./255; % Color
    c_4 = [2, 53, 53]./255; % Color
    
    patch(x_1, y, c_1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_2, y, c_2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_3, y, c_3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    patch(x_4, y, c_4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Add Text!
    x_center_1 = mean(x_1(1:2));
    x_center_2 = mean(x_2(1:2));
    x_center_3 = mean(x_3(1:2));
    x_center_4 = mean(x_4(1:2));
    y_center = mean(y);
    
    text(x_center_1, y_center, '$\mathcal{S}_{1}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_2, y_center, '$\mathcal{S}_{2}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_3, y_center, '$\mathcal{S}_{3}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    text(x_center_4, y_center, '$\mathcal{S}_{4}$', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle')
    
    xline(2, '-.')
    xline(4, '-.')
    xline(6, '-.')
    
    % Power_Comparison(1,:) = t;
    % Power_Comparison(2,:) = P_I_1;
    % Power_Comparison(3,:) = P_I_2;
    % Power_Comparison(4,:) = P_I_3;
    
    p4 = plot(Temp_Violation(1,:), Temp_Violation(2,:), 'Color', ...
        '#d95555', 'LineWidth', 1.5);
    
    p5 = plot(Temp_Violation(1,:), Temp_Violation(3,:), 'Color', ...
        '#7EB283', 'LineWidth', 1.5);
    
    p6 = plot(Temp_Violation(1,:), Temp_Violation(4,:), 'Color', ...
        '#1a46a3', 'LineWidth', 1.5);
    hold off
    
    T_max_1 = max(Temp_Violation(2,:));
    T_max_2 = max(Temp_Violation(3,:));
    T_max_3 = max(Temp_Violation(4,:));
    T_max = max([T_max_1, T_max_2, T_max_3]);
    
    xlim([0, 8])
    ylim([0, 8])
    
    xlabel('Time in h')
    ylabel('$\hat{l}_\vartheta$ in $^{\circ}C$s')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    
    x_section_bar = [1, 2, 3, 4];
    
    b_1 = bar(x_section_bar, P_vals_bar);
    b_1(1).FaceColor = '#d95555';
    b_1(2).FaceColor = '#7EB283';
    b_1(3).FaceColor = '#1a46a3';
    
    
    % legend('$C_1$','$C_2$', '$C_3$')
    xticks([1 2 3 4])
    xticklabels({'$\mathcal{S}_{1}$','$\mathcal{S}_{2}$', ...
        '$\mathcal{S}_{3}$', '$\mathcal{S}_{4}$'})
    xlabel('Section')
    y2 = ylabel('$\hat{P}$ in $A^2$s');
    % y2.Position(1) = -0.4;
    label_pos_2 = get(y2, 'Position');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    
    b_2 = bar(x_section_bar, T_vals_bar);
    b_2(1).FaceColor = '#d95555';
    b_2(2).FaceColor = '#7EB283';
    b_2(3).FaceColor = '#1a46a3';
    
    xticks([1 2 3 4])
    xticklabels({'$\mathcal{S}_{1}$','$\mathcal{S}_{2}$', ...
        '$\mathcal{S}_{3}$', '$\mathcal{S}_{4}$'})
    xlabel('Section')
    ylabel('$\hat{l}_\vartheta$ in $^{\circ}C$s')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    
    b_3 = bar(P_vals_bar, 'stacked');
    b_3(1).FaceColor = [15, 194, 192]./255;
    b_3(2).FaceColor = [0, 143, 140]./255;
    b_3(3).FaceColor = [1, 89, 88]./255;
    b_3(4).FaceColor = [2, 53, 53]./255;
    
    xticks([1 2 3])
    xticklabels({'$\mathcal{C}_{1}$', ...
        '$\mathcal{C}_{2}$','$\mathcal{C}_{3}$'})
    xlabel({'Controller', 'a) Energy consumption'})
    y3 = ylabel('$\hat{P}$ in $A^2$s');
    y3.Position(1) = -0.65;
    label_pos_3 = get(y3, 'Position');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    
    b_4 = bar(T_vals_bar, 'stacked');
    b_4(1).FaceColor = [15, 194, 192]./255;
    b_4(2).FaceColor = [0, 143, 140]./255;
    b_4(3).FaceColor = [1, 89, 88]./255;
    b_4(4).FaceColor = [2, 53, 53]./255;
    
    xticks([1 2 3])
    xticklabels({'$\mathcal{C}_{1}$', ...
        '$\mathcal{C}_{2}$','$\mathcal{C}_{3}$'})
    xlabel({'Controller', 'b) Temperature violation'})
    ylabel('$\hat{l}_\vartheta$ in $^{\circ}C$s')
    
    
    leg = legend([b_3(1), b_3(2), b_3(3), b_3(4), b_2(1), b_2(2), b_2(3)], ...
        {'Section $\mathcal{S}_{1}$','Section $\mathcal{S}_{2}$','Section $\mathcal{S}_{3}$', ...
        'Section $\mathcal{S}_{4}$', 'Controller $\mathcal{C}_{1}$', 'Controller $\mathcal{C}_{2}$', ...
        'Controller $\mathcal{C}_{3}$', 'Door openings'}, 'NumColumns', 2);
    
    leg.Title.String = '\rm Bar charts';
    
    set(leg,'visible','on');
    set(leg, 'box', 'off');
    set(leg, 'Orientation', 'Vertical');
    leg.Layout.Tile = 'north';
    
    leg2 = legend([p1, p2, p3, r], ...
        {'Controller $\mathcal{C}_{1}$', 'Controller $\mathcal{C}_{2}$', ...
        'Controller $\mathcal{C}_{3}$', 'Door opening'}, 'NumColumns', 1);
    
    leg2.Title.String = '\rm Line charts';
    
    set(leg2,'visible','on');
    set(leg2, 'box', 'off');
    set(leg2, 'Orientation', 'Vertical');
    set(leg2, 'Interpreter', 'latex')
    leg2.Layout.Tile = 'north';
    
    set_figure_properties(hfig_efficiency)

    if save_efficiency == true
        saveas(gcf,'Efficiency','pdf');
    end
end



%% Q_z Analysis of Cargo load

if show_Qz == true

    load('configuration_z.mat');
    load('configuration_z0.mat');
    load('Q_z_4_smooth.mat');

    mode_1 = [255, 224, 140]./255; % #FFE08C
    mode_3 = [255, 163, 140]./255; % #FFA38C darker: #D99977
    c_z = [0, 128, 128]./255;

    t_z = configuration_z.tout;
    T_ICB_z = configuration_z.T_ICB.Data;
    s_z = configuration_z.s_z.Data;
    T_cargo = configuration_z.T_cargo.Data;
    I_z = configuration_z.I.Data;

    T_ICB_z0 = configuration_z0.T_ICB.Data;
    I_z0 = configuration_z0.I.Data;

    % Cut from 4th phase:
    idx_z_4 = find(t_z == 6*3600);
    t_z_4 = t_z(idx_z_4:end);

    t_z_4 = t_z_4 - t_z_4(1);
    t_z_4 = t_z_4/3600;

    T_ICB_z4 = T_ICB_z(idx_z_4:end);
    s_z4 = s_z(idx_z_4:end);
    I_z4 = I_z(idx_z_4:end);

    T_ICB_z0 = T_ICB_z0(idx_z_4:end);
    I_z0 = I_z0(idx_z_4:end);


    hfig_z = figure;

    tiledlayout(9, 1,'TileSpacing','compact','Padding','tight');

    nexttile([2, 1])

    door_controller_z = door_controller(:,end);
    door_controller_z = door_controller_z - 6;

    x = [door_controller_z(1,end), door_controller_z(2,end), door_controller_z(2,end), ...
        door_controller_z(1,end)];
    x_1 = [0, door_controller_z(1,end), door_controller_z(1,end), 0];
    x_3 = [door_controller_z(2,end), 6, 6 door_controller_z(2,end)];
    y = [0, 0,  1000, 1000];
    c = door_color;
    m2 = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on
    m1 = patch(x_1, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    m3 = patch(x_3, y, mode_3, 'FaceAlpha', 0.5, 'EdgeColor', 'none');


    p1 = plot(t_z_4, Q_z_4_smooth, 'Color', c_z, 'LineWidth', 1.5);
    hold off

    xticks([0, 1, 2, 3, 4])
    xlim([0, 4])
    ylim([0 800])

    y1 = ylabel('$\dot{Q}_{z}$ in W');
    % y_1.Position(1) = -0.4;
    label_pos_1 = get(y1, 'Position');

    nexttile([2, 1])

    for j = 1:numel(T_ICB_z4)
        if T_ICB_z4(j) >= 15
            T_ICB_z4(j) = 15;
        end
    end


    x_hyst = [0, door_controller_z(1,end), door_controller_z(1,end), 0];
    y_hyst = [4.5, 4.5, 5.5, 5.5];
    c_hyst = [.5 .5 .5]; % Color
    hyst = patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    hold on

    x_hyst = [door_controller_z(2,end), 6, 6, door_controller_z(2,end)];
    patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, 'EdgeColor', 'none')

    x = [door_controller_z(1,end), door_controller_z(2,end), door_controller_z(2,end), ...
        door_controller_z(1,end)];
    y = [0, 0,  1000, 1000];
    c = door_color;
    patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    y = [5.5, 5.5,  1000, 1000];
    patch(x_1, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    patch(x_3, y, mode_3, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    y = [0, 0,  4.5, 4.5];
    patch(x_1, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    patch(x_3, y, mode_3, 'FaceAlpha', 0.5, 'EdgeColor', 'none');


%     x = [door_controller_z(1,end), door_controller_z(2,end), door_controller_z(2,end), ...
%         door_controller_z(1,end)];
%     y = [0, 0,  1000, 1000];
%     c = door_color;
%     r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%     hold on
%     m1 = patch(x_1, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%     m3 = patch(x_3, y, mode_3, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    p2 = plot(t_z_4, T_ICB_z4, 'Color', '#A60A33',  'LineWidth', 1.5);
    hold off
  

    xticks([0, 1, 2, 3, 4])
    xlim([0, 4])
    ylim([4, 15])

    y2 = ylabel('$\vartheta_{cc}$ in $^{\circ}C$');
    y2.Position(1) = -0.21;
    label_pos_2 = get(y2, 'Position');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nexttile([2, 1])

    for i = 1:numel(T_ICB_z0)
        if T_ICB_z0(i) >= 15
            T_ICB_z0(i) = 15;
        end
    end


    x_hyst = [0, door_controller_z(1,end), door_controller_z(1,end), 0];
    y_hyst = [4.5, 4.5, 5.5, 5.5];
    c_hyst = [.5 .5 .5]; % Color
    hyst = patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    hold on

    x_hyst = [door_controller_z(2,end), 6, 6, door_controller_z(2,end)];
    patch(x_hyst, y_hyst, c_hyst, 'FaceAlpha', 0.2, 'EdgeColor', 'none')

    x = [door_controller_z(1,end), door_controller_z(2,end), door_controller_z(2,end), ...
        door_controller_z(1,end)];
    y = [0, 0,  1000, 1000];
    c = door_color;
    patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    y = [5.5, 5.5,  1000, 1000];
    patch(x_1, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    patch(x_3, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    y = [0, 0,  4.5, 4.5];
    patch(x_1, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    patch(x_3, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    p2 = plot(t_z_4, T_ICB_z0, 'Color', '#A60A33',  'LineWidth', 1.5);

    dummy_1 = yline(-1000, 'Color', [1, 1, 1]);
    dummy_2 = yline(-2000, 'Color', [1, 1, 1]);
    hold off

    xticks([0, 1, 2, 3, 4])
    xlim([0, 4])
    ylim([4, 15])

    y2 = ylabel('$\vartheta_{cc}$ in $^{\circ}C$');
    y2.Position(1) = -0.21;
    label_pos_2 = get(y2, 'Position');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nexttile([2, 1])

    x = [door_controller_z(1,end), door_controller_z(2,end), door_controller_z(2,end), ...
        door_controller_z(1,end)];
    y = [0, 0,  1000, 1000];
    c = door_color;
    r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on
    m1 = patch(x_1, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    m3 = patch(x_3, y, mode_3, 'FaceAlpha', 0.5, 'EdgeColor', 'none');


    p3 = plot(t_z_4, I_z4, 'LineWidth', 1.5);
    hold off

    xticks([0, 1, 2, 3, 4])
    xlim([0, 4])
    ylim([0, 7])

    y3 = ylabel('$I_{tec}$ in A');
    y3.Position(1) = -0.21;
    label_pos_3 = get(y3, 'Position');

    nexttile

    x = [door_controller_z(1,end), door_controller_z(2,end), door_controller_z(2,end), ...
        door_controller_z(1,end)];
    y = [0, 0,  1000, 1000];
    c = door_color;
    r = patch(x, y, c, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on
    m1 = patch(x_1, y, mode_1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    m3 = patch(x_3, y, mode_3, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    p4 = plot(t_z_4, s_z4, 'Color', c_z, 'LineWidth', 1.5);
    hold off

    xlim([0, 4])
    ylim([0, 1])
    xticks([0, 1, 2, 3, 4])
    yticks([0 1])
    y4 = ylabel('$s_z$ in 1');
    y4.Position(1) = -0.21;
    label_pos_4 = get(y4, 'Position');

    xlabel('Time in h')

    leg = legend([m1, m2, m3, hyst, dummy_1, dummy_2, p2, p1, p3], ...
        {'Mode 1','Mode 2','Mode 3', 'Hysteresis band', '', '', ...
        'System outputs', ...
        'Disturbances', 'System inputs'}, 'NumColumns', 3);

    set(leg,'visible','on');
    set(leg, 'box', 'off');
    set(leg, 'Orientation', 'Vertical');
    leg.Layout.Tile = 'north';

    set_figure_properties(hfig_z)

    if store_Qz == true
        saveas(gcf,'Q_z','pdf');
    end
end



%% Picture Functions:

function set_figure_properties(hfig)
    picturewidth = 21; % 21 for A4 width!
    hw_ratio = 0.7; % Figure Height (changeable) 1.8!!!
    % Note: 1.3 or 0.85 or 0.5
    
    % Set figure properties
    set(findall(hfig, '-property','FontSize'),'FontSize',14)
    set(findall(hfig, '-property','Box'),'Box','off') % optional (off)
    set(findall(hfig, '-property','Interpreter'),'Interpreter','latex')
    set(findall(hfig, '-property','TickLabelInterpreter'), ...
        'TickLabelInterpreter','latex')
    set(hfig, 'Units','centimeters','Position', ...
        [0 0 picturewidth hw_ratio*picturewidth])
    pos = get(hfig, 'Position');
    set(hfig, 'PaperPositionMode','Auto', 'PaperUnits','centimeters', ...
        'PaperSize',[picturewidth, hw_ratio*picturewidth]) % manual - Auto
%     print(hfig, 'pdf_figure', '-dpdf', '-vector', '-fillpage');


end



function hOut = TextLocation(textString,varargin,dataset)
    
    t = annotation('textbox');
    t.String = textString;
    t.LineStyle = 'None';
    
    % Set properties to match the rest of the plot
    t.FontSize = 12;
    t.Interpreter = 'latex';
    if dataset == "val"
        t.Color = [102,205,170]./255;
    else
        t.Color = [153 153 255]./255;
    end
    % Check if a valid position vector was provided
    if numel(varargin) == 1 && isnumeric(varargin{1}) && numel(varargin{1}) == 4
        t.Position = varargin{1};
    else
        % Use default position if no valid position vector was provided
        t.Position = varargin;
%         t.Position = [0.5 0.5 0 0];
    end
    
    if nargout
        hOut = t;
    end
end


% Uncomment to get the location!
% function hOut = TextLocation(textString,varargin)
%     
%     l = legend(textString,varargin{:});
%     t = annotation('textbox');
%     t.String = textString;
%     t.Position = l.Position
%     set(l,'visible','off');
% %     delete(l);
%     t.LineStyle = 'None';
%     
%     % Set properties to match the rest of the plot
%     t.FontSize = 17;
%     t.Interpreter = 'latex';
%     
%     if nargout
%         hOut = t;
%     end
% end


