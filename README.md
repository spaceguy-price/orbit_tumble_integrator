# Orbit Generation 
## Overview
Project imports orbital information and generates rigid body tumble motion using a 4th Order Runge-Kutta Integrator. The code outputs two files, the tumble attitude and the sun vector. Both are expressed in quaternion.


|<img src=../../docs/media/orbit_LLA.png width=300>|<img src=../../docs/media/orbit_ECI.png width=300>|
|---|---|
|<img src=../../docs/media/sunV_Hill.png width=300>|<img src=../../docs/media/omega.png width=300>|

# Requirements
|Matlab Version| R2023a|
|---|---|
|Toolbox | Aerospace |
||Signal Processing|

## Operation Procedure

To run the code, open Q_sunV_generator.m and modify the settings section found at line 38.

Uses one function from the MATLAB Aerospace toolbox to convert between reference frames.
### Orbit Setting
```matlab
% Orbital Settings
ID = 33500;                             % [int] TLE object identifier. See ./utils/getSatelliteTLE.m to add different orbits.
dtUTC = datetime('now', 'TimeZone','Z');% dtUTC = datetime(2023, 5, 25, 12, 0, 0, 'TimeZone', 'Z'). The date influences position in orbit and sun vector calculation.
```
* Select predefined TLE represented orbit by __ID__
* Or, you can define orbit as TLE [getSatelliteTLE.m](./utils/getSatelliteTLE.m)
* __dtUTC__ is a time setting parameter to decide sun light vector and on-orbit position


### Sun Vector Setting
```matlab
number_orbits = 1;                      % [double]  Number of orbits to simulate
sunV_ECI = approxECISunPosition(dtUTC); % Sun Vector Settings [x,y,z] unit vector. NOTE: KHI project special request sun vector can be found in section "Obtain sun vector" (~line 200)
sunV_ECI = sunV_ECI/norm(sunV_ECI);     % Or enter your own vector. Example: sunV = [1,0,0];
```
* In the code above, __sunV_ECI__, which is sun vector in Earth Centered Inertia Frame, is automatically calculated by __dtUTC__.
* Or, you can set sun vector in ECI directory. They are expressed in [x, y, z] elements. 
### Tumbling Setting
```matlab
% Tumble Motion Settings
Q_import = false;                       % [boolean] If importing tumble motion data set to true
% If importing tumble data (Q_import = true)
import_filename = 'gravity_gradient_motion_array.mat'; % File path of the data to be imported
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
```
|Time[s]|In Plane Angle[deg]|Out Of Plane Angle[deg]|Pitch[rad]|Roll[rad]|Target Nose Vector in Hill Frame x[m]|y[m]|z[m]|
|---|---|---|---|---|---|---|---|
|0|45|-45|0.785|0.615|-0.707|0.408|0.577|
|...|...|...|...|...|...|...|...|
|1999|-30|-44|-0.529|0.696|0.504|0.553|0.663|

* You can chose import motion data or simulate in this code. 
* If you set __Q_import__ as true, .mat file will be imported. This is array which looks like above. (not included headers)
* If you set __Q_import__ as false, you need set initial state of target. 
* __q_initial__ is initial attitude quaternion. Default is [1, 0, 0, 0] and it is rarely modified. 
* __w_initial__ is initial velocity of target. First element should be zero, and others represent rotation around [x, y, z]
* __I__ is initial inertia matrix 
