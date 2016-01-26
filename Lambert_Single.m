%% Robert Miller
clc
clear;
format compact

% Information:
% Developed by Robert Miller - UAF Student
% Uses the lambert process
% Solves
    % 2 Velocites (v1, v2)
% Required
    % Functions - NONE
    % r1
    % r2
    % delta_T0 = Time of Flight
    
    
%% Given Information: 
r1 = [0.1385, 0.4187, 0.8948];          % DU
r2 = [1.025097, 1.653562, 0.801700];    % DU
delta_T0 = 10/13.446849;                % TU    (Time of Flgiht)
mu = 1;                                 % DU^3/TU^2
tol = 10^-15;                           % Tolerance
v2_S = [-0.143047, 0.243260, -0.169915];   % DU/TU


fprintf('Given Information\n\n')
fprintf('r1 (DU) \t\t = \t [%.4f, %.4f, %.4f] DU\n', r1(1), r1(2), r1(2))
fprintf('r2 (DU) \t\t = \t [%.4f, %.4f, %.4f] DU\n', r2(1), r2(2), r2(3))
fprintf('Delta T0 (TU)\t = \t %.4f TU\n',delta_T0)


% Notes:
% 1 DU = 6378 km
% 1 TU = 13.446849 s

%% Calculations:

r1_norm = norm(r1);
r2_norm = norm(r2);

% Solving or cos(anomaly) and sin(anomaly)
cos_delta_anomaly = dot(r1,r2)/(r1_norm*r2_norm);
sin_delta_anomlay = abs(cross(r1,r2))/(r1_norm*r2_norm);


% Solvin A (Note: if A = 0, then unable to solve problem)
% For short way
tm = 1;

% For long way
% tm = -1

A = tm*sqrt(r1_norm*r2_norm*(1+cos_delta_anomaly));     % Correct, based on exmaple
fprintf('A \t\t\t\t = \t %.4f\n', A)

% Used to diplay what type of orbit
if (A == 0)
   fprintf('\n\nERROR: A = 0\nUnable to compute orbital elements')
elseif(A ~= 0)
    if(tm == 1)
        fprintf('\nOrbital Calculation: Short Way (tm = +1)\n')
    elseif(tm == -1)
       fprintf('\nOrbital Calculation: Long Way (tm = -1)\n')
    end
end


% Counter
counter = 1;

% Initial Starting Conditions:
z = 0;
C = 1/2;
S = 1/6;
z_up = 4*pi^2;
z_low = -4*pi;


fprintf('\n\n------------------------------------------------------\n')
fprintf('Iteration\t\ty\t\t\tx\t\t\tDelta_T\t\t\tz\n\n')
while (counter <= 10000)
    
    % Prof. [C,S] Function
    if (z >= 1e-6)
        C = (1-cos(z^0.5))/z;
        S = (z^.5-sin(z^0.5))/z^1.5;
    elseif (z <= -1e-6)
        C = (1-cosh((-z)^0.5))/z;
        S = (sinh((-z)^0.5)-(-z)^0.5)/(-z)^1.5;
    else
        C = 0.5;
        S = 1/6;
    end
    
    y = r1_norm + r2_norm + A*(z*S-1)/sqrt(C);      % Correct for first iteration
    x = sqrt(y/C);                                  % Correct for first iteration
    delta_T = ((x^3)*S + A*sqrt(y))/sqrt(mu);       % Correct for first iteration

    
    if(delta_T <= delta_T0)
        z_low = z;
    else
        z_up = z;
    end
    
    fprintf('%.f', counter)
    fprintf('\t\t\t\t%.4f', y)
    fprintf('\t\t%.4f', x)
    fprintf('\t\t%.4f', delta_T)
    fprintf('\t\t\t%.4f\n', z)
    
    % Breaks Loop
    ratio = abs(delta_T - delta_T0);
    if(ratio <= tol)
        break;
    end
    
    % Sets new value of z based on z_up and z_low
    z = (z_up + z_low)/2;
        
    %z = z+1            % Appears to be causing a bug when adding +1
    counter = counter+1;
    
    % Debug Tool
    % fprintf('\n\n==========|  Next Iteration  |==========\n\n')
end

%fprintf('\n===============================================================\n')
fprintf('\nIteration Completed: \nUnits in (AU)\n')
fprintf('y \t\t\t = \t %.8f\n', y)
fprintf('x \t\t\t = \t %.8f\n', x)
fprintf('z \t\t\t = \t %.8f\n', z)
fprintf('Delta_T \t = \t %.8f TU\n', delta_T)
fprintf('Delta_T0 \t = \t %.8f TU\n', delta_T0)
fprintf('Error \t\t = \t %.8f percent\n\n', ratio*100)

fprintf('\n===============================================================\n')
fprintf('Calculating Universal Variables: \n')

% Calculating Universal Variables
f = 1 - (y/r1_norm);
g = A*sqrt(y/mu);
g_dot = 1 - (y/r2_norm);

v1_I = (r2-f*r1)/(g);         % Units (DU/TU)
v2_I = (g_dot*r2 - r1)/(g);   % Units (DU/TU)

fprintf('f \t\t\t = \t %.8f\n', f)
fprintf('g \t\t\t = \t %.8f\n', g)
fprintf('g_dot \t\t = \t %.8f\n\n', g_dot)

fprintf('v1 (DU/TU) \t = [%.5f, %.5f, %.5f] DU/TU\n', v1_I(1), v1_I(2), v1_I(3))
fprintf('v1 (km/s) \t = [%.5f, %.5f, %.5f] km/s\n\n', v1_I(1)*7.90536828, ...
    v1_I(2)*7.90536828, v1_I(3)*7.90536828)

fprintf('v2 (DU/TU) \t = [%.5f, %.5f, %.5f] DU/TU\n', v2_I(1), v2_I(2), v2_I(3))
fprintf('v2 (km/s) \t = [%.5f, %.5f, %.5f] km/s\n\n', v2_I(1)*7.90536828, ...
    v2_I(2)*7.90536828, v2_I(3)*7.90536828)


% Post Calculations:
delta_V1 = (v2_I - v1_I);
delta_V2 = (v2_S - v2_I);
delta_V = delta_V1 + delta_V2;
fprintf('delta_V1 (Intercepter)\t= [%.4f, %.4f, %.4f]\t DU/TU\n',delta_V1(1), delta_V1(2), delta_V1(3))
fprintf('delta_V2 (Intercepter) \t= [%.4f, %.4f, %.4f]\t DU/TU\n',delta_V2(1), delta_V2(2), delta_V2(3))
fprintf('delta_V (Total) \t\t= [%.4f, %.4f, %.4f]\t DU/TU\n\n\n',delta_V(1), delta_V(2), delta_V(3))


fprintf('\n===============================================================\n')
delta_V1a = norm(delta_V1);
delta_V2a = norm(v2_S - v2_I);
delta_V = delta_V1a+delta_V2a;

fprintf('delta_V1 (Intercepter)\t= %.4f\t DU/TU\n',delta_V1a)
fprintf('delta_V2 (Intercepter)\t= %.4f\t DU/TU\n',delta_V2a)
fprintf('delta_V2 (Intercepter)\t= %.4f\t DU/TU\n',delta_V)









