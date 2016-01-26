%% Robert W. Miller

clear
clc
format compact
format shortG

%% Notes

        % Calculating position coordinates in IJK (relative to earth) with
        % earth being considered ellipsoidal and not spherical.
         
        
        % ae = Equatorial Radius (km) - use if object is in this vicinity
        % be = Polar Radius (km) - use if object is in this vicinity
        
        % L = Latitude of ground station (degrees)
        % H = Height/Elevation of object (km)
        % e = Eccentricity of object
        % theta_LST = Local Side Real Time (degrees)
        
        
%% Given Information

% Note
    % 1 DU = 6378.1 km
    % 1 DU/TU = 7.90538 km/s
    % 1 TU = 806.80415 sec

% Earth Radius
ae = 6378.145;      % km  |  Equatorial Radius 
be = 6356.785;      % km  |  Polar Radius


% Input Information
L = 64.84;
H = 720;     
e = .08182;
theta_LST = 195.2883;
%% Calculations

x_sat = ((ae*cosd(L)) / sqrt(1 - e^2*sind(L)^2) ) + H*cosd(L);

z_sat = ((ae*(1-e^2)*sind(L) ) / (sqrt(1-e^2*sind(L)))) + H*sind(L);
r_IJK_km = [x_sat*cosd(theta_LST), x_sat*sind(theta_LST), z_sat];
r_IJK_DU = [x_sat*cosd(theta_LST), x_sat*sind(theta_LST), z_sat]/6378.145;

%% Printing Calculations
fprintf('---------- INPUT INFORMTION ----------\n')
fprintf('L \t\t\t= %.4f \t\tdegrees\n', L);
fprintf('H \t\t\t= %.4f \t\tkm\n', H);
fprintf('e \t\t\t= %.4f\n', e);
fprintf('Theta LST \t= %.4f \t\tdegrees\n\n\n', theta_LST);

fprintf('---------- CALCULATED INFORMATION ----------\n')
fprintf('Position IJK = [%.4f, %.4f, %.4f] km\n', r_IJK_km(1), r_IJK_km(2), r_IJK_km(3))
fprintf('Position IJK = [%.4f, %.4f, %.4f] DU\n\n', r_IJK_DU(1), r_IJK_DU(2), r_IJK_DU(3))



