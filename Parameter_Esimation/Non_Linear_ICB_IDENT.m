function [dxdt, y] = Non_Linear_ICB_IDENT(t,x,u,theta_1,theta_2,theta_3,...
    theta_4,theta_5,theta_6,theta_7,theta_8,theta_9,theta_10,theta_11,...
    theta_12,theta_13,varargin)
% varargin ... immer notwendig als Zusatz für die Durchführung der
%              Identifikation

global c % bekannte Systemparameter aus globaler Variable laden
global f_1 % For the first Iteration
%global f_2 % For the second Iteration

% unbekannte Systemparameter in Vektorform schreiben - Reihenfolge
% beachten!
theta = [theta_1;theta_2;theta_3;theta_4;theta_5;theta_6;theta_7;...
    theta_8;theta_9;theta_10;theta_11;theta_12;theta_13]; 
% Before: theta_13;theta_14  

[dxdt, ~] = Non_Linear_ICB(t,x,c,f_1,theta,u,[]);



% Festlegung des Ausgangsvektor, d.h. was gemessen wurde und für die
% Identifikation als Systemausgang zur Verfügung steht.
y = [x(3);...
     x(4)];

end
