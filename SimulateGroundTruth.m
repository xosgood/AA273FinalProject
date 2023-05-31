% 5/22/23
% AA273 Final Project
% Simulate our 2D leader bird dynamics and 
clc; clear; close all;

rng(273);

n_L = 3; % number of dimensions of state of leader
n_F = 3; % number of dimensions of state of follower
p = 2; % number of dimensions of control input

dt = 0.1;
t_f = 200;
tspan = 0:dt:t_f;
N = length(tspan);

Q = 0.001 * eye(n_L);
Q_rel = 0.2 * eye(n_F);

v = 10 * ones(1,N); % velocity command
omega = sin(tspan/3); % angular velocity command
u_L = [v; omega];

u_F = zeros(p, N);
v_follower_max_thresh = 5;
omega_follower_max_thresh = 0.05;

K_p = [.1; .01] * 10 * sin(5*tspan); % proportional gain for control loop for follower

x_L = zeros(n_L, N); % true leader state, through simulating dynamics
x_F_des = zeros(n_F, N); % desired follower state, through simulating dynamics
x_F_act = zeros(n_F, N); % actual follower state

x_L(:,1) = zeros(n_L, 1);
x_F_des(:,1) = [-10; -10; deg2rad(0)];
x_F_act(:,1) = [-20; -20; deg2rad(45)];

%% simulate dynamics
for i = 2:N
    w = mvnrnd(zeros(n_L,1), Q)'; % process noise
    w_rel = mvnrnd(zeros(n_F,1), Q_rel)'; % process noise
    x_L(:,i) = f_abs(x_L(:,i-1), u_L(:,i-1), dt) + w; % leader dynamics propagation
    
    x_F_des(:,i) = f_rel(x_F_des(:,i-1), [0; 0], dt); % desired follower dynamics propagation
    x_F_act(:,i) = f_rel(x_F_act(:,i-1), u_F(:,i-1), dt) + w_rel; % actual follower dynamics propagation
    
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
end

x_F_des_global = x_L + x_F_des;
x_F_act_global = x_L + x_F_act;


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
figure; grid on; hold on;
plot(x_L(1,:), x_L(2,:))
plot(x_F_des_global(1,:), x_F_des_global(2,:), '.');
plot(x_F_act_global(1,:), x_F_act_global(2,:), 'o');
xlabel("x"); ylabel("y");
title("Leader bird trajectory");
legend("Leader","Follower desired","Follower actual");


%% functions
% nonlinear leader dynamics
function x_new = f_abs(x_old, u, dt)
    x_new = zeros(size(x_old));
    x_new(1) = x_old(1) + dt * u(1) * cos(x_old(3));
    x_new(2) = x_old(2) + dt * u(1) * sin(x_old(3));
    x_new(3) = x_old(3) + dt * u(2);
end

% nonlinear measurement
function y = g(x)
    y = norm(x(1:2));
end

% desired relative dynamics
function x_new = f_rel(x_old, u, dt)
    x_new = zeros(size(x_old));
    x_new(1) = x_old(1) + dt * u(1) * cos(x_old(3));
    x_new(2) = x_old(2) + dt * u(1) * sin(x_old(3));
    x_new(3) = x_old(3) + dt * u(2);
end
