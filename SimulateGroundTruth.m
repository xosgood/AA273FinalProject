% 5/22/23
% AA273 Final Project
% Simulate our 2D leader bird dynamics and 
clc; clear; close all;

rng(273);

n_L = 3; % number of dimensions of state of leader
p = 2; % number of dimensions of control input

dt = 0.1;
t_f = 100;
tspan = 0:dt:t_f;
N = length(tspan);

Q = 0.001 * eye(n_L);

v = 10 * ones(1,N); % velocity command
omega = sin(tspan/3); % angular velocity command
u = [v; omega];

x = zeros(n_L, N); % true state, through simulating dynamics

mu_0 = 0; % initial state estimate
Sigma_0 = 0.01 * eye(n_L); % initial covariance estimate

%% simulate dynamics
for i = 2:N
    w = mvnrnd(zeros(n_L,1), Q)'; % process noise
    x(:,i) = f(x(:,i-1), u(:,i-1), dt) + w; % dynamics propagation
end

% plot
figure;
plot(tspan, x(1,:), tspan, x(2,:), tspan, x(3,:))
figure;
plot(x(1,:), x(2,:))


%% functions
% nonlinear dynamics
function x_new = f(x_old, u, dt)
    x_new = zeros(size(x_old));
    x_new(1) = x_old(1) + dt * u(1) * cos(x_old(3));
    x_new(2) = x_old(2) + dt * u(1) * sin(x_old(3));
    x_new(3) = x_old(3) + dt * u(2);
end

% nonlinear measurement
function y = g(x)
    y = norm(x(1:2));
end

