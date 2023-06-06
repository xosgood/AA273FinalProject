% 5/22/23
% AA273 Final Project
% Simulate our 2D leader bird dynamics and 
clc; clear; close all;

rng(273);
%State Definitions
n_L = 3; % number of dimensions of state of leader
n_F = 3; % number of dimensions of state of follower

%Control and Measurement Definition
p = 2; % number of dimensions of control input
m = 3; % number of dimensions of measurement (leader's measurement of follower)

%Run Time
t_f = 1; %[sec]
dt = 0.1;

tspan = 0:dt:t_f;
N = length(tspan);

%Probability for Process Noise 
Q = 0.001 * eye(n_L);    % For Leader
Q_rel = 0.2 * eye(n_F);  % For Follower

%Probability for Measurement Noise
R = diag([1, 0.5, 0.5]); %For Leader 

%
v = 10 * ones(1,N);   % velocity command
omega = sin(tspan/3); % angular velocity command
u_L = [v; omega];

u_F = zeros(p, N);
v_follower_max_thresh = 5; % Max Velocity Change Per Time Step
omega_follower_max_thresh = 0.05; %Max Angle Change Per Time Step

K_p = [.1; .01] * 10 * sin(5*tspan); % proportional gain for control loop for follower

% creating vectors to hold data
x_L = zeros(n_L, N); % true leader state, through simulating dynamics
x_F_des = zeros(n_F, N); % desired follower state, through simulating dynamics
x_F_act = zeros(n_F, N); % actual follower state
y = zeros(m, N-1); % measurements


% initial conditions
x_L(:,1) = zeros(n_L, 1);
x_F_des(:,1) = [-10; 10; deg2rad(0)]; %Bird
x_F_act(:,1) = [-20; 20; deg2rad(45)];


follower_obj = Follower(n_F, t_f, dt, [-10; -10; deg2rad(0)], [-12; -12; deg2rad(5)], Q_rel);
follower_obj1 = Follower(n_F, t_f, dt, [-10; 10; deg2rad(0)], [-12; 12; deg2rad(-5)], Q_rel);

w_rel_arr = mvnrnd(zeros(follower_obj.n_F,1), follower_obj.Qtrue, follower_obj.numsteps-1)';

%% simulate dynamics
for i = 2:N
    % leader dynamics
    w = mvnrnd(zeros(n_L,1), Q)'; % process noise
    x_L(:,i) = f_abs(x_L(:,i-1), u_L(:,i-1), dt) + w; % leader dynamics propagation
    
    % follower dynamics
  
    follower_obj.curr_ind = i;
    follower_obj.desiredDynamics(); 
    follower_obj.actualDynamics(w_rel_arr(:, i-1)); 
    follower_obj.generateControl();
    
    
    follower_obj1.curr_ind = i;
    follower_obj1.desiredDynamics(); 
    follower_obj1.actualDynamics(w_rel_arr(:, i-1)); 
    follower_obj1.generateControl();
    
    w_rel = w_rel_arr(:, i-1);%mvnrnd(zeros(n_F,1), Q_rel)'; % process noise
    x_F_des(:,i) = f_rel(x_F_des(:,i-1), [0; 0], dt); % desired follower dynamics propagation
    x_F_act(:,i) = f_rel(x_F_act(:,i-1), u_F(:,i-1), dt) + w_rel; % actual follower dynamics propagation
    
    
    % control law for follower
    e = x_F_act(:,i-1) - x_F_des(:,i-1);
    e_x = e(1);
    e_y = e(2);
    e_psi = e(3);
    e_rho = norm(e(1:2));
    Beta = atan2(e_y, e_x);
    u_F(:,i) = K_p(:,i) .* [e_rho;
                            e_psi - Beta];
    u_F(1,i) = sign(u_F(1,i)) * min(abs(u_F(1,i)), v_follower_max_thresh);
    u_F(2,i) = sign(u_F(2,i)) * min(abs(u_F(2,i)), omega_follower_max_thresh);
    
    u_obj = follower_obj.u_F(:, i);
    u_loc = u_F(:, i);   
end

x_F_des_global = x_L + x_F_des;
x_F_act_global = x_L + x_F_act;

x_F_des_global1 = x_L + follower_obj.x_F_des;
x_F_act_global1 = x_L + follower_obj.x_F_act;

x_F_des_global2 = x_L + follower_obj1.x_F_des;
x_F_act_global2 = x_L + follower_obj1.x_F_act;



%% plotting
% plot leader
figure; grid on; hold on;
plot(tspan, x_L(1,:), tspan, x_L(2,:), tspan, x_L(3,:))
xlabel("time (s)"); ylabel("pose");
title("Leader bird pose over time");
legend("x","y","\theta");
figure; grid on; hold on;
plot(x_L(1,:), x_L(2,:));
xlabel("x"); ylabel("y");
title("Leader bird trajectory");

% plot follower
figure; grid on; grid minor; hold on;
plot(x_L(1,:), x_L(2,:), 'b-', 'LineWidth', 2)
plot(x_F_des_global1(1,:), x_F_des_global1(2,:), 'r.');
plot(x_F_act_global1(1,:), x_F_act_global1(2,:), 'o', 'MarkerSize', 10, 'color', [255, 165, 0]/255); %'MarkerFaceColor', [255, 165, 0]/255
plot(x_F_des_global2(1,:), x_F_des_global2(2,:), 'g.');
plot(x_F_act_global2(1,:), x_F_act_global2(2,:), 'mo', 'MarkerSize', 10); %'MarkerFaceColor', 'm'
xlabel("x"); ylabel("y");
view([90 -90])
title("Three Ship Formation Flying");
legend("Leader", "Desired Right Side Follower", "Right Actual Follower", "Desired Left Side Follower", "Left Actual Follower");

figure; 
title("Bird Trajectories");
subplot(2,1,1)
hold on; grid on; grid minor;
plot(tspan, x_F_des_global(1,:), 'b');
plot(tspan, x_F_act_global(1,:), 'r');
%plot(tspan, x_F_act_global1(1, :), 'g');
ylabel("X Position");
hold off;

subplot(2,1,2);
hold on; grid on; grid minor;
plot(tspan, x_F_des_global(2,:), 'b');
plot(tspan, x_F_act_global(2,:), 'r');
plot(tspan, x_F_act_global1(2, :), 'g');
ylabel("Y Position"); 
xlabel("Time [s]");
legend("Follower Desired", "Follower Actual", "OOP Actual", "Location", "South");
hold off;

% figure;
% subplot(2,1,1)
% plot(tspan, u_F(1,:), 'b', 'LineWidth', 1.5);
% hold on; grid on; grid minor;
% plot(tspan, follower_obj.u_F(1,:), 'r');
% ylabel("Velocity Command");
% hold off;
% 
% subplot(2,1,2)
% plot(tspan, u_F(2,:), 'b', 'LineWidth', 1.5);
% hold on; grid on; grid minor;
% plot(tspan, follower_obj.u_F(2,:), 'r');
% ylabel("Theta Command");
% legend("Control", "OOP Cntrl", "Location", "Best");
% hold off;


%% functions
% nonlinear leader dynamics
function x_new = f_abs(x_old, u, dt)
    x_new = zeros(size(x_old));
    x_new(1) = x_old(1) + dt * u(1) * cos(x_old(3));
    x_new(2) = x_old(2) + dt * u(1) * sin(x_old(3));
    x_new(3) = x_old(3) + dt * u(2);
end

% nonlinear measurement
% function y = g(x)
%   
% end

% desired relative dynamics
function x_new = f_rel(x_old, u, dt)
    x_new = zeros(size(x_old));
    x_new(1) = x_old(1) + dt * u(1) * cos(x_old(3));
    x_new(2) = x_old(2) + dt * u(1) * sin(x_old(3));
    x_new(3) = x_old(3) + dt * u(2);
end
