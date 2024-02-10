function [dxdt, add_out] = Non_Linear_ICB(t,x,c,f_1,theta,u,t_vec)

% ------------------- ERKLÄRUNG DER STRUKTUR -----------------------------
% Outputs:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% dxdt      ... zeitliche Ableitung des Zustandsvektors
% add_out   ... zusätzliche auszugebende Größen aus der Funktion
%
% Inputs:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% t         ... Zeit 
% x         ... Zustandsvektor
% u         ... Eingangsvektor
% d         ... Vektor nicht beeinflussbarer, aber messbarer Eingänge
% z         ... Störungsvektor
% c         ... Vektor der konstanten, bekannten Systemparameter
% theta     ... Vektor der konstanten, unbekannten Systemparameter
% ------------------------------------------------------------------------

% ----------------- ERKLÄRUNG DER MODELLGRÖSZEN --------------------------
% Zustände:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        

% x(1) = T_wall_inside          Temperature [K]
% x(2) = T_wall_outside         Temperature [K]      
% x(3) = T_cargo                Temperature [K]
% x(4) = T_hs                   Temperature [K]
% x(5) = T_ICB                  Temperature [K]
%
% Konstante, bekannte Systemparameter:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% c(1) = First Parameter for the outer wall temperature             
% c(2) = Second Parameter for the outer wall temperature         
%
% Konstante, unbekannte Systemparameter:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% theta(1) = a1        Temperaturkonstante a1 in 1/s
% theta(2) = a2        Temperaturkonstante a1 in 1/s
% theta(3) = a3        Temperaturkonstante a1 in 1/s
% theta(4) = a4        Temperaturkonstante a1 in 1/s
% theta(5) = a6        Temperaturkonstante a1 in 1/s
% theta(6) = a7        Temperaturkonstante a1 in V/J
% theta(7) = b2        Temperaturkonstante a1 in (V*K)/(A*J)
% theta(8) = a8        Temperaturkonstante a1 in 1/s
% theta(9) = a9        Temperaturkonstante a1 in 1/s
% theta(10) = a11      Temperaturkonstante a1 in 1/s
% theta(11) = xi       Temperaturkonstante a1 in 1/s
% theta(12) = b3       Temperaturkonstante a1 in K/s
%
% Eingangsgrößen:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% u(1) = I              Strom am Peltier Element in A
% u(2) = s_DOOR          Logical Vector for Fan
% u(3) = s_FAN         Logical Vector for Door Opening
% u(4) = T_Water (Fresh Water In) Temperature in Kelvin
% u(5) = T_amb (ambient Temperature) Temperature in Kelvin
%
% messbare, nicht veränderliche Eingangsgrößen:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% d = []                hier NICHT vorhanden!
%
% Störgrößen:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% z(1) = T_water (u(4))                  Störgrösse Umgebungstemperatur in K
% z(2) = T_amb (u(5))                    Störgrösse Wassertemperatur in K
%
% ------------------------------------------------------------------------

if ~isempty(t_vec)
    u = (interp1(t_vec,u',t,'previous'))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_18 = par.Data;
% t_Data = Data_18(:,1);
% 
% s_Door_ref = Data_18(:,6);
% s_Door = interp1(t_Data, s_Door_ref, t, 'previous');
% 
% I_ref = Data_18(:,8);
% I = interp1(t_Data, I_ref, t, 'previous');
% 
% s_Fan_ref = Data_18(:,7);
% s_Fan = interp1(t_Data, s_Fan_ref, t, 'previous');

% Before adding u(4) & u(5) this was the original implementation, but it
% is not needed any more:
% % T_atm = c(1);
% % T_water = c(2);

%Beschreibung des Eingangvektors u:
% u(1) = I              Strom am Peltier Element in A
% u(2) = s_DOOR          Logical Vector for Fan
% u(3) = s_FAN         Logical Vector for Door Opening
% u(4) = T_Water (Fresh Water In) Temperature in Kelvin
% u(5) = T_amb (ambient Temperature) Temperature in Kelvin

%Beschreibung des Vektors x:
% x(1) -> T_wall_inside
% x(2) -> T_wall_outside
% x(3) -> T_cargo
% x(4) -> T_hs
% x(5) -> T_ICB

% NOTE:
% First Iteration only theta(i) * f_1(i)
% Second Iteration: theta(i) * f_1(i) * f_2(i)

%All 5 Differential Equations
dxdt(1) = theta(1)*f_1(1)*(x(2)-x(1)) - ...
    theta(2)*f_1(2)*(x(1)-x(4));
dxdt(2) = theta(3)*f_1(3)*(u(5)-x(2)) - c(2)*(x(2)-x(1)); % All parameters fixed!
% dxdt(3) = c(1)*(x(3)-x(5)); % Equation for Cargo
dxdt(3) = (theta(4)*f_1(4)*u(3) + ...
    (1-u(3))*theta(5)*f_1(5))* ...
    (x(4)-x(3)) - ...
    theta(6)*f_1(6)*u(1)*x(3) +...
    theta(7)*f_1(7)*(u(1)^2) + ...
    theta(8)*f_1(8)*(u(4)-x(3)); 
dxdt(4) = theta(9)*f_1(9)*(x(1)-x(4)) - ...
    (theta(10)*f_1(10)*u(3) + ...
    (1-u(3))*theta(11)*f_1(11))* ... % PROBLEM HERE!!!! (Klammer setzen)
    (x(4)-x(3)) + ...
    theta(12)*f_1(12)*(u(5)-x(4))*u(2) + ...
    theta(13)*f_1(13)*u(3); 
    % + c(2)*(x(3)-x(5)) % If Cargo needs to be implemented

dxdt = dxdt';


% Festlegung zusätzlicher Funktionsausgänge, die für gewisse
% Problemstellungen interessant sind
% Simulink kann mit [] als Ausgang nicht arbeiten, wodurch in diesem Fall
% der Ausgang zu 0 zu setzen wäre.
add_out = [x(3);...
           x(4)];

end
