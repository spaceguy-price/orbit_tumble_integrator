%% Q_sunV_generator.m
% Code to generate tumble motion quaternions and sun unit vectors for an
% object in a predefined orbit about Earth. This code has been heavily
% modified from original code designed to import TLE data and showcase the
% orbit in p-q, ECI and LLA formats.

%  When       Who    What
%  ---------- ------ --------------------------------------------------
%  2019/07/11 mnoah  original code (orbit and different frames)
%  2021/12/14 aprice included tumble integrator and sun-angle progression
%  2023/02/02 aprice updated for phase 2 of the KHI project
%  2023/05/22 aprice commented code for handoff to SRL orbital

% As of 10 Dec, 2021 the sun is 147,329,779 km from Earth
% https://theskylive.com/how-far-is-sun
% GOSAT approximate altitude = 677 km
% GOSAT orbital inclination at 98.1 degrees
% GOSAT is in a sun-synchronous orbit
% GOSAT period is 98.1 minutes

% Taking sun location (sub-solar point) from https://rl.se/sub-solar-point
% Can also obtain atronomy information from https://rhodesmill.org/pyephem/
%
% pip install pyephem
% python
% >>> import ephem
% >>> s = ephem.Sun()
% >>> s.earth_distance

% Include dependency subfolders containing additional functions
addpath('./utils');
addpath('./tumble_integrator');
close all % close old figures
% The standard gravitational parameter $$\mu$$ of a celestial body is the
% product of the gravitational constant G and the mass M of the body. 
mu = 3.98618e14; % [m3/s2] Earth's geocentric gravitational constant

%% Settings
%--------------------------------------------------------------------------
%------------------------- SETTINGS START ---------------------------------
%--------------------------------------------------------------------------
% Update this section before running the code.

% Orbital Settings
ID = 33500;                             % [int] TLE object identifier. See ./utils/getSatelliteTLE.m to add different orbits.
dtUTC = datetime('now', 'TimeZone','Z');% dtUTC = datetime(2023, 5, 25, 12, 0, 0, 'TimeZone', 'Z'). The date influences position in orbit and sun vector calculation.
number_orbits = 1;                      % [double]  Number of orbits to simulate
sunV_ECI = approxECISunPosition(dtUTC); % Sun Vector Settings [x,y,z] unit vector. NOTE: KHI project special request sun vector can be found in section "Obtain sun vector" (~line 200)
sunV_ECI = sunV_ECI/norm(sunV_ECI);     % Or enter your own vector. Example: sunV = [1,0,0];
orbit_plot = true;                      % If you want MATLAB to generate three different plots visualizing the orbit

% Tumble Motion Settings
Q_import = false;                       % [boolean] If importing tumble motion data set to true
% If importing tumble data (Q_import = true)
import_filename = 'KHI_gravity_gradient_motion_array.mat'; % File path of the data to be imported
import_time_step = 1;                   % [s] The imported data has a timestep that may be different than our requested timestep
% If generating tumble data (Q_import = false)
q_initial = [1,0,0,0]';                 % Initial attitude quaternion [qw, qx, qy, qz]' 
w_initial = [0,0,0,pi()/180]';          % [0, wx, wy, wz]' [rad/s] initial relative rotation
dt = 0.1;                               % [s] simulation timestep (4th order Runge-Kutta) <--- NOTE: time_step >= dt
sim_plot = true;                        % [boolean] whether to show the energy and angular velocity plots
sim_viewer = false;                     % [boolean] whether to show the orientation viewer (NOTE: the program will not continue while the viewer is running)
I = [27548.8, 0,       0;               % [kgm2] object inertia matrix
     0,       27675.3, 0; 
     0,       0,       6808.8]; 

clear Ixx Iyy Izz Ixy Iyz Ixz

% File Settings
time_step = 10;                         % [s] output timestep. NOTE: The tumble integrator will output data with interval dt. We then downsample to time_step. Be careful when importing tumble data to match timesteps.
Q_filename = 'targetQ_10s.txt';         % quaternion savefile path
sunV_filename = 'sunV_10s.txt';         % sunvector savefile path
%--------------------------------------------------------------------------
%------------------------- SETTINGS END -----------------------------------
%--------------------------------------------------------------------------
%% Obtain Orbital Elements
[TLE] = getSatelliteTLE(ID);

[OE] = TLE2OrbitalElements(TLE);
fprintf(1,['Kepler Elements for satelliteID %d epoch %s:\n' ...
    '\ta [m] = %f\n\te = %f\n\ti [deg] = %f\n\tomega [deg] = %f\n' ...
    '\tOmega [deg] = %f\n\tM [deg] = %f\n'], floor(OE.satelliteID), ...
    datestr(OE.epoch), OE.a_km*1e3, OE.e, OE.i_deg, OE.omega_deg, ...
    OE.Omega_deg, OE.M_deg);
a_m = OE.a_km*1e3;
M_deg = OE.M_deg;

%% *Orbital Plane Coordinate Frame*
%  p_m - [m] coordinate along axis through center and perigee
%  q_m - [m] coordinate passing through focus and perpendicular to p-axis
%  dpdt_m_per_s = [rad/s] p component velocity
%  dqdt_m_per_s = [rad/s] q component velocity
n_rad_per_s = sqrt(mu/a_m^3);  % [rad/s] mean motion
n_deg_per_s = rad2deg(n_rad_per_s); % [deg/s] mean motion
M_rad = deg2rad(M_deg);
E_rad = M_rad; 
dE = 99999;
eps = 1e-6; % [rad] control precision of Newton's method solution
while (abs(dE) > eps)
    dE = (E_rad - OE.e * sin(E_rad) - M_rad)/(1 - OE.e * cos(E_rad));
    E_rad = E_rad -  dE;
end
p_m = a_m*(cos(E_rad) - OE.e);
q_m = a_m*sqrt(1 - OE.e^2)*sin(E_rad);

% Orbital angular velocity variables in p-q frame (aprice doesn't use)
% dMdt_rad_per_s = n_rad_per_s;
% dEdt_rad_per_s = dMdt_rad_per_s/(1 - e*cos(E_rad));
% h = -a_m*sin(E_rad)*dEdt_rad_per_s;
% dqdt_m_per_s = a_m*cos(E_rad)*dEdt_rad_per_s*sqrt(1 - e^2);
E_deg_epoch = rad2deg(E_rad); 

n_orbits = str2double(OE.n_orbits_per_day);
OE.orbital_period = 24*60*60/n_orbits;          %[s] Length of a single orbit

clear n_orbits dEdt_rad_per_s dpdt_m_per_s dMdt_rad_per_s h 

%% Break-up orbit into time steps and rotate orbital data To ECI
Rz_Omega = [ ...
    [cosd(OE.Omega_deg) sind(OE.Omega_deg) 0]; ...
    [-sind(OE.Omega_deg) cosd(OE.Omega_deg) 0]; ...
    [0 0 1]];
Rx_i = [ ...
    [1 0 0]; ...
    [0 cosd(OE.i_deg) sind(OE.i_deg)]; ...
    [0 -sind(OE.i_deg) cosd(OE.i_deg)]];
Rz_omega = [ ...
    [cosd(OE.omega_deg) sind(OE.omega_deg) 0]; ...
    [-sind(OE.omega_deg) cosd(OE.omega_deg) 0]; ...
    [0 0 1]];

% time of epoch
[Year,Month,Day,H,M,S] = datevec(OE.epoch);
date = datetime(Year, Month, Day, H, M, S);
jd = juliandate(date);
T = (jd - 2451545.0d0)./36525.0d0; % form time in Julian centuries from J2000
Theta_deg = 100.460618375 + 36000.770053608336*T + 0.0003879333*T^2 + 15*H + M/4 + mod(S/240,360);

Rz_hour = [ ...
    [cosd(Theta_deg) sind(Theta_deg) 0]; ...
    [-sind(Theta_deg) cosd(Theta_deg) 0]; ...
    [0 0 1]];

% position of the satellite at epoch in the orbit p-q plane
r_pq = [p_m q_m 0]';
% position of satellite at epoch in ECI coordinates
omega_deg = OE.omega_deg;
Omega_deg = OE.Omega_deg;
i_deg = OE.i_deg;
r_ECI = inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pq;
%r_LLA = eci2lla(r_ECI',datevec(datenum(Year, Month, Day, H, M, S)),'IAU-2000/2006'); % Requires matlab AEROSPACE toolbox
[r_ECEF, ~, ~] = ECItoECEF(jd , r_ECI, [0;0;0], [0;0;0]);
[r_LLA] = ecef2lla(r_ECEF');

da = (360/OE.orbital_period)*time_step;     % [deg] orbital angular progression for "time_step" defined in settings
Evals = 0:da:360*number_orbits;             % increment through the orbit(s) in "time_step"[s] intervals
Orbit_p = a_m*(cosd(Evals)-OE.e);              % [m] orbit positions p
Orbit_q = a_m*sqrt(1 - OE.e^2)*sind(Evals);    % [m] orbit positions q
deltaT_s = ((Evals-E_deg_epoch) - OE.e*sind(Evals-E_deg_epoch))/n_deg_per_s; % [s] time since epoch along orbit

% Convert p-q frame orbit data to ECI and LLA frames.
% NOTE: The ECItoECEF and ecef2lla functions are vectorized and thus we could forgo the for loop
Orbit_ECI = zeros(numel(deltaT_s),3);
Orbit_LLA = zeros(numel(deltaT_s),3);
for ipt = 1:size(Orbit_ECI,1)
    r_pq = [Orbit_p(ipt) Orbit_q(ipt) 0]';
    Orbit_ECI(ipt,:) = (inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pq)';
    [temp, ~, ~] = ECItoECEF(jd , Orbit_ECI(ipt,:)', [0;0;0], [0;0;0]);
    [lla] = ecef2lla(temp');
    Orbit_LLA(ipt,:) = lla;
end

clear ipt Evals Year Month Day H M S da Orbit_p Orbit_q deltaT_s dt0 date jd
clear a_m dE E_rad

%% Obtain Hill Frame in ECI coordinates for each orbital position
% x-z define the orbital plane with -z pointing to the Earth's center
% The y axis is therefore coincident with the orbital plane normal vector
X_orbitF = zeros(size(Orbit_ECI));
Z_orbitF = X_orbitF;
eul = [OE.Omega_deg, 0, OE.i_deg,]*(pi()/180);
Y_orbitF = ones(size(Orbit_ECI)).* (eul2rotm(eul,'ZYX')*[0;0;1])';

for ipt = 1:size(Orbit_ECI,1)
    Z_orbitF(ipt,:) = Orbit_ECI(ipt,:)./(norm(Orbit_ECI(ipt,:)));
    X_orbitF(ipt,:) = cross(Y_orbitF(ipt,:), Z_orbitF(ipt,:));
end

X_hillF = Z_orbitF;             % Convert custom orbital frame to the Hill frame
Y_hillF = X_orbitF;
Z_hillF = Y_orbitF;

clear eul ipt X_orbitF Y_orbitF Z_orbitF

%% Obtain sun vector in the Hill frame
% Create special case sun vectors (KHI Project Request)
%{
% Z_hillF is the orbital plane normal vector
Z_ = Z_hillF(1,:);
% Rotate Y_ECI about Z_ECI to get a specific angle
% Solving for the rest is a solving a quadratic equation
% 45 degrees, we can arbitrarily set sun_ECIz = 0
input45 = 45;
p(1) =  (1/Z_(1))^2 * (Z_(2)^2) + 1;
p(2) = -2 * (1/Z_(1))^2 * Z_(2) * cosd(input45);
p(3) =  (1/Z_(1))^2 * (cosd(input45))^2 - 1;
r = roots(p);
sun45(2) = r(1);                % pick one of the roots for the y component
sun45(3) = 0;                   % set the z component to zero
sun45(1) = sqrt(1-sun45(2)^2);  % solve for the x component
sun45 = sun45./norm(sun45);     % convert to a unit vector
sun45_angle = 90 - atan2d(norm(cross(Z_,sun45)), dot(Z_,sun45));    %check that the angle the sun vector makes with the orbital plane is correct
%}

% Map the sun vector to each Hill frame
sunV = zeros(length(X_hillF),3);
for ipt = 1:size(sunV,1)
    R = [X_hillF(ipt,:); Y_hillF(ipt,:); Z_hillF(ipt,:)];
    
    sunV(ipt,:) = (R*(sunV_ECI'))';
    sunV(ipt,:) = sunV(ipt,:)./norm(sunV(ipt,:));
    
end

clear p r ipt R n_orbits input45 Z_

%% Generate or Import Target Motion

% Generate Motion
if Q_import == false
    sim_t = number_orbits*OE.orbital_period;    %[s] simulation duration
    % NOTE: the sim_viewer coordinate frame does not necessarily match the frame defined by KHI
    [targetQ] = tumble_integrator(q_initial, w_initial, I, dt, sim_t, sim_plot, sim_viewer);

% Import Motion
else
    % 1) Load the file
    load(import_filename)

    % 2) Separate the data into vectors
    t = KHI_data(:,1);          %[s] Time
    theta_in = KHI_data(:,2);   %[deg] In-plane angle
    theta_out = KHI_data(:,3);  %[deg] Out-of-plane angle
    pitch = KHI_data(:,4);      %[rad] Pitch angle
    roll = KHI_data(:,5);       %[rad] Roll Angle
    v_x = KHI_data(:,6);        %[m] Target nose vector in cartesian coordinates (Hill Frame)
    v_y = KHI_data(:,7);
    v_z = KHI_data(:,8);
    ze = zeros(length(v_x),1);

    v_ = [v_x, v_y, v_z];

    % Convert KHI frame to Blender Hill Frame
    % Rotate 90 degs positive about the y axis
    % Rotate 180 degs about the x axis
    RotY = [cos(-pi/2),0,sin(-pi/2); 0,1,0; -sin(-pi/2),0,cos(-pi/2)]; 
    RotX = [1,0,0; 0,cos(pi),-sin(pi); 0,sin(pi),cos(pi)];
    v_ = (RotX'*RotY'*v_')';

    % 3) Convert [cartesian] to [spherical]
    [inp, outp, r] = cart2sph(v_(:,1), v_(:,2), v_(:,3));
    % check that r remains constant
    check = (r - 5.3) > 0.001;
    if any(check)
        disp('The attitude input data was not interpreted correctly.')
    end

    % 4) Create attitude quaternions based on the yaw and pitch angles
    Q = eul2quat([inp, -outp, ze], 'ZYX');         % Good but adds a slight visual roll
    %Q = eul2quat([yaw, ze, pi/2 - pitch], 'XYZ'); % Contains a LOT of roll

    % 5) Shorten the file to be the same length as three orbits
    targetQ = quaternion(Q(1:uint32(OE.orbital_period/import_time_step),:));    %<--- NOTE: This depends on <number_orbits>, <time_step> and the KHI <import_time_step>
end

clear r check RotY RotZ v_x v_y v_z KHI_data

%% Write to file
% Downsample the tumble motion to match the requested timestep
if Q_import
    if import_time_step > time_step
       error('The requested time step is smaller than the imported tumble motion time step.')
    end
    ds_targetQ = downsample(targetQ(1:end), time_step/import_time_step);
    
else
    if dt > time_step
        error('The requested time step is smaller than the tumble motion integrator timestep.')
    end
    ds_targetQ = downsample(targetQ(1:end), time_step/dt);
end

% Write attitude to file
% [qw, qx, qy, qz]
[qA,qB,qC,qD] = parts(ds_targetQ(1:end));
fid=fopen(Q_filename,'wt');
fprintf(fid, '%f %f %f %f \n', [qA qB qC qD]');
fclose(fid); 

% Write sun unit vectors to file
% [x, y, z]
fid=fopen(sunV_filename,'wt');
fprintf(fid, '%f %f %f \n', sunV');
fclose(fid);

clear fid qA qB qC qD ds_targetQ

%% Plots
if orbit_plot

e_r = 6378137.0; % Earth's approximate equatorial radius
    
% Plot the sun vector as it appears in the Hill frames
figure
plot3([0,40*sunV(1,1)], [0,40*sunV(1,2)], [0,40*sunV(1,3)], 'k-', 'LineWidth', 4)
hold on
grid on
plot3(0,-40,0,'ko')
for i = 2:uint32(3*length(sunV)/4)
    plot3([0,40*sunV(i,1)], [0,40*sunV(i,2)], [0,40*sunV(i,3)], 'b-')
end
xlabel('Hill X [m]');
ylabel('Hill Y [m]');
zlabel('Hill Z [m]');
title('Sun Vector Check (Hill Frame)');
axis equal

% *Plot Cartesian Coordinates* (Earth-Centered Inertial)
figure
hold on

% Orbit
plot3(Orbit_ECI(:,1), Orbit_ECI(:,2), Orbit_ECI(:,3), 'm-', 'linewidth', 4)% orbit
plot3(Orbit_ECI(1,1), Orbit_ECI(1,2), Orbit_ECI(1,3), ...                  % starting point of orbit
    'p','MarkerFaceColor','k','MarkerEdgeColor','none','Markersize',30)

% Sun Vector
plot3([sunV_ECI(1)*e_r*3, 0], [sunV_ECI(2)*e_r*3, 0], [sunV_ECI(3)*e_r*3, 0], ':', 'linewidth', 10, 'color', [0.8500 1*0.25 0.0980]);

% Hill frame axes
ax_sz = 0.5;        % length of Hill frame axes
for i = uint32(linspace(1,length(Orbit_ECI),5))
    plot3([Orbit_ECI(i,1), Orbit_ECI(i,1) + X_hillF(i,1)*e_r*ax_sz], ...
          [Orbit_ECI(i,2), Orbit_ECI(i,2) + X_hillF(i,2)*e_r*ax_sz], ...
          [Orbit_ECI(i,3), Orbit_ECI(i,3) + X_hillF(i,3)*e_r*ax_sz], ...
          'r-', 'linewidth', 3);
    plot3([Orbit_ECI(i,1), Orbit_ECI(i,1) + Y_hillF(i,1)*e_r*ax_sz], ...
          [Orbit_ECI(i,2), Orbit_ECI(i,2) + Y_hillF(i,2)*e_r*ax_sz], ...
          [Orbit_ECI(i,3), Orbit_ECI(i,3) + Y_hillF(i,3)*e_r*ax_sz], ...
          'g-', 'linewidth', 3);
    plot3([Orbit_ECI(i,1), Orbit_ECI(i,1) + Z_hillF(i,1)*e_r*ax_sz], ...
          [Orbit_ECI(i,2), Orbit_ECI(i,2) + Z_hillF(i,2)*e_r*ax_sz], ...
          [Orbit_ECI(i,3), Orbit_ECI(i,3) + Z_hillF(i,3)*e_r*ax_sz], ...
          'b-', 'linewidth', 3);
end
% Earth object
[x,y,z] = sphere;
x = x*e_r;
y = y*e_r;
z = z*e_r;
surf(x,y,z)

% Labels
xlabel('ECI x [m]');
ylabel('ECI y [m]');
zlabel('ECI z [m]');
title('Target in Earth Centered Inertial (ECI) Coordinates');
grid on
axis equal
legend('Target Orbit', 'Hill Frame Axes', 'SunV', 'Hill X', 'Hill Y', 'Hill Z');

% Plot orbit path on map of the Earth
figure('color','white');
earth = imread('ear0xuu2.jpg');
lv= size(earth,1);
lh= size(earth,2);
lats =  (1:lv)*180/lv - 90;
lons =  (1:lh)*360/lh - 180;
image(lons, -lats, earth)
hold on;
set(gca,'ydir','normal');
grid on
plot(Orbit_LLA(:,2),Orbit_LLA(:,1),'.r', 'Markersize', 12);
plot(r_LLA(2),r_LLA(1),'p','MarkerFaceColor',[1 0 0.7],'MarkerEdgeColor','none','Markersize',30);
set(gca,'XTick',[-180:30:180]);
set(gca,'YTick',[-90:30:90]);
set(gca,'Xcolor',0.3*ones(1,3));
set(gca,'Ycolor',0.3*ones(1,3));
title(['Target Ground Track for Epoch ' datestr(OE.epoch)])
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
axis equal
axis tight
set(gca,'fontweight','bold');

end