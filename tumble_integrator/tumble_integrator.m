function [Q] = tumble_integrator(q, w, I, dt, sim_t, plot_on, view_on)
% Rigid Body Dynamics 4th Order Rune-Kutta Integrator
% Andrew Price
% 2021 December 15

% Inputs:
%   q - initial attitude quaterionion [q0,qx,qy,qz] 
%   w - initial rotation vector (body-inertial) [rad/s] [0,x,y,z]
%   I - object inertia matrix 3x3 [kgm2]
%   dt - simulation timestep [s]. Recommend a small timestep for stability (ex 0.01s)
%   sim_t - simulation duration [s]
%   view_on - [boolean], a 1 will launch the viewer 

%% Simulation Parameters
% System total energy (used to check if the integration is valid)
T_total = 0.5*dot(w(2:end), I*w(2:end));

% Simulation parameters
t = 0:dt:sim_t;
N = length(t);

% Preallocate Output
T = zeros(N,1);
W = zeros(N,3);
Q = quaternion(T,T,T,T);

%% 4th Order Runge Kutta Integration Simulation (RK4)
    
% Euler's euqations of motion
f = @(ww, I) -inv(I)*cross(ww, I*ww);

% Initialize k vectors (slope estimators)
k1 = zeros(3,1);
k2 = zeros(3,1);
k3 = zeros(3,1);
k4 = zeros(3,1);

for i = 1:N
    % (1) Update attitude
    qR = quatDE(q);
    q_dot = 0.5*qR*w;
    q = q + q_dot*dt;

    q = q/sum(q.^2);    % Enforce quaternion constraint: sum(qi^2) = 1
    Q(i) = quaternion(q');

    % (2) Runge Kutta Update Acceleration
    w_ = w(2:end);  %rotation rate x,y,z [rad/s]

    k1 = f(w_            , I);
    k2 = f(w_ + (dt/2)*k1, I);
    k3 = f(w_ + (dt/2)*k2, I);
    k4 = f(w_ + dt*k3    , I);

    w_ = w_ + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % (3) Calculate energy state
    T(i) = 0.5*dot(w_, I*w_);

    % (4) Update rotation vector
    w(2:end) = w_;
    W(i,:) = w_';

    % Display progress
    if mod(t(i),100) == 0
        str1 = strcat('Progress:', int2str(t(i)), '/', int2str(sim_t), ' [s] ');
        str2 = strcat(' || Energy: ', sprintf('%.3f', T(i)*1000), ' [mJ]');
        disp(strcat(str1, str2))
    end

end

%% Visualize
%Assumptions
%-Estimated inertias
%-Inertia Z axis aligned with spacecraft body axisymmetric axis
%-No external torques or gravity gradient

if plot_on == 1
    figure
    t = [0:dt:sim_t]';
    plot(t(1:end-100), W(1:end-100,1), 'r-');
    hold on
    set(gca,'fontsize',20)
    title(strcat('Progression of \omega', 'FontSize', 30));
    xlabel('time [s]', 'FontSize', 24);
    ylabel('Rotation Rate [rad/s]', 'FontSize', 24);
    plot(t, W(:,2), 'g-');
    plot(t, W(:,3), 'b-');
    legend('\omega_X', '\omega_Y', '\omega_Z');
    hold off

    figure
    plot(t(1:end-100),1000*T(1:end-100))
    title('Energy Progression Versus Time', 'FontSize', 30)
    xlabel('Time [s]', 'FontSize', 24)
    ylabel('Energy [mJ]', 'FontSize', 24)
end

if view_on == 1
    viewer_x = 100;   % Set the speed of the replay
    string = strcat('Tumbling Motion [', int2str(viewer_x),'X replay]');
    viewer = Orientation_Viewer('Title',{string});
    for i= 1:length(Q)
        viewer(Q(i))
        pause(dt/viewer_x)
    end
end