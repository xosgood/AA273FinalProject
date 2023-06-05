% 5/22/23
% AA273 Final Project
% Simulate our 2D leader bird dynamics and 
clc; clear; close all;

rng(273);

a = 1.96; % for 95% confidence interval

n_L = 3; % number of dimensions of state of leader
num_followers = 2; % number of follower birds
n_F = 3 * num_followers; % number of dimensions of state of follower
p = 2 * num_followers; % number of dimensions of control input
m = 3 * num_followers; % number of dimensions of measurement (leader's measurement of follower)

dt = 0.1;
t_f = 100;
tspan = 0:dt:t_f;
N = length(tspan);

Q_abs = 0.001 * eye(n_L);
Q = 1 * eye(n_F);
R = 2 * diag(repmat([0.1, 0.01, 0.01], 1, num_followers));

v = 10 * ones(1,N); % velocity command
omega = sin(tspan/3); % angular velocity command
u_L = [v; omega];

u_F = zeros(p, N);
v_follower_thresh = 5;
omega_follower_thresh = 0.05;
u_F_thresh = [v_follower_thresh; omega_follower_thresh];

% proportional gain for control loop for follower
%K_p = [.1; .01] * 10 * sin(5*tspan);
K_cycle = 1000;
K_p1 = repmat([0.001; 0.01], num_followers, 1); % good
K_p2 = repmat([0.1; 0.01], num_followers, 1); % bad


% creating vectors to hold data
x_L = zeros(n_L, N); % true leader state, through simulating dynamics
x_F_des = zeros(n_F, N); % desired follower state, through simulating dynamics
x_F_act = zeros(n_F, N); % actual follower state
y = zeros(m, N-1); % measurements
mu = zeros(n_F, N); % state estimate
Sigma = zeros(n_F, n_F, N); % covariance estimate

% initial conditions
x_L(:,1) = zeros(n_L, 1);
x_F_des(:,1) = [-10; -15; deg2rad(0); -10; 15; deg2rad(0)];
x_F_act(:,1) = [-20; -25; deg2rad(45); -20; 25; deg2rad(-45)];
mu(:,1) = -ones(n_F, 1);
Sigma(:,:,1) = eye(n_F);

% timing
loop_times = zeros(1, N-1);

%% loop
for i = 2:N
    %%% Ground truth
    % leader dynamics
    w_abs = mvnrnd(zeros(n_L,1), Q_abs)'; % process noise
    x_L(:,i) = f(x_L(:,i-1), u_L(:,i-1), dt) + w_abs; % leader dynamics propagation
    
    % follower dynamics
    w = mvnrnd(zeros(n_F,1), Q)'; % process noise
    x_F_des(:,i) = f(x_F_des(:,i-1), zeros(p, 1), dt); % desired follower dynamics propagation
    x_F_act(:,i) = f(x_F_act(:,i-1), u_F(:,i-1), dt) + w; % actual follower dynamics propagation
    
    e_for_controller = zeros(p, 1);
    p_iter = 1;
    for j = 1:3:n_F
        % control law for follower (control for next timestep)
        e = x_F_act(j:j+2,i-1) - x_F_des(j:j+2,i-1);
        e_x = e(1);
        e_y = e(2);
        e_psi = e(3);
        e_rho = norm(e(1:2));
        Beta = atan2(e_y, e_x);
        e_for_controller(p_iter:p_iter+1) = [e_rho;
                                   e_psi - Beta];
        p_iter = p_iter + 2;
    end
    if mod(i, 2 * K_cycle) <= K_cycle
        K_p = K_p1;
    else
        K_p = K_p2;
    end
    u_F(:,i) = K_p .* e_for_controller;
    u_F(:,i) = sign(u_F(:,i)) .* min(abs(u_F(:,i)), repmat(u_F_thresh, num_followers, 1));
    
    %%% Measure (leader's measurments of followers)
    v = mvnrnd(zeros(m,1), R)'; % measurement noise
    y(:,i-1) = g(x_F_act(:,i)) + v; % measure
    
    %%% EKF
    t1 = tic;
    
    A = DynamicsJacobian(mu(:,i-1), u_F(:,i-1), dt); % Jacobian for dynamics
    mu(:,i) = f(mu(:,i-1), u_F(:,i-1), dt); % mean predict
    Sigma(:,:,i) = A * Sigma(:,:,i-1) * A' + Q; % covariance predict
    C = MeasurementJacobian(mu(:,i)); % Jacobian for measurements
    K = Sigma(:,:,i) * C' / (C * Sigma(:,:,i) * C' + R); % Kalman gain
    mu(:,i) = mu(:,i) + K * (y(:,i-1) - g(mu(:,i))); % mean update
    Sigma(:,:,i) = Sigma(:,:,i) - K * C * Sigma(:,:,i); % covariance update
    
    loop_times(i-1) = toc(t1);
end

avg_loop_time = mean(loop_times)

x_F_des_global = repmat(x_L, num_followers, 1) + x_F_des;
x_F_act_global = repmat(x_L, num_followers, 1) + x_F_act;

% find confidence intervals around relative state
tspan_conf = [tspan tspan(end:-1:1)];
mu_conf = [mu(1,:) + a * sqrt(reshape(Sigma(1,1,:), 1,length(Sigma(1,1,:)))), mu(1,end:-1:1) - a * sqrt(reshape(Sigma(1,1,end:-1:1), 1,length(Sigma(1,1,end:-1:1))));
           mu(2,:) + a * sqrt(reshape(Sigma(2,2,:), 1,length(Sigma(2,2,:)))), mu(2,end:-1:1) - a * sqrt(reshape(Sigma(2,2,end:-1:1), 1,length(Sigma(2,2,end:-1:1))));
           mu(3,:) + a * sqrt(reshape(Sigma(3,3,:), 1,length(Sigma(3,3,:)))), mu(3,end:-1:1) - a * sqrt(reshape(Sigma(3,3,end:-1:1), 1,length(Sigma(3,3,end:-1:1))))];


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

% plot follower and measurements
figure; grid on; hold on;
plot(tspan, x_F_des(1,:), tspan, x_F_des(2,:));
plot(tspan, x_F_act(1,:), tspan, x_F_act(2,:));
plot(tspan(2:end), y(1,:) .* y(2,:), "o", tspan(2:end), y(1,:) .* y(3,:), "o");
xlabel("time (s)"); ylabel("pose");
title("Follower bird true relative state over time");
legend("Desired x", "Desired y", "Actual x", "Actual y", "Measured x", "Measured y");

% plot follower true state and estimated state
figure;
sgtitle("Follower bird true & estimated relative state over time \n with 95% confidence interval shaded in");
subplot(2,1,1); grid on; hold on;
plot(tspan, x_F_act(1,:), tspan, x_F_act(2,:));
plot(tspan, mu(1,:), tspan, mu(2,:));
p = fill(tspan_conf, mu_conf(1,:),'yellow');
%p.FaceColor = [0.8 0.8 1];
p.FaceAlpha = 0.25;
p.EdgeColor = 'none';
p = fill(tspan_conf, mu_conf(2,:),'m');
%p.FaceColor = [1 0.8 0.8];
p.FaceAlpha = 0.25;
p.EdgeColor = 'none';
sgtitle("Follower bird true and estimated relative position over time");
xlabel("time (s)"); ylabel("pose");
legend("True x", "True y", "Estimated x", "Estimated y", "95% confidence interval on estimated x", "95% confidence interval on estimated y");
subplot(2,1,2); grid on; hold on;
plot(tspan, x_F_act(3,:));
plot(tspan, mu(3,:));
p = fill(tspan_conf, mu_conf(3,:),'r');
%p.FaceColor = [0.8 0.8 1];
p.FaceAlpha = 0.25;
p.EdgeColor = 'none';
xlabel("time (s)"); ylabel("relative heading [rad]");
title("Follower bird true and estimated relative heading over time");
legend("True \psi", "Estimated \psi", "95% confidence interval on estimated \psi");


%% functions
% nonlinear dynamics (treating as non-holonomic robot)
function x_new = f(x_old, u, dt)
    x_new = zeros(size(x_old));
    j = 1;
    for i = 1:3:length(x_old)
        v = u(j); % velocity command
        w = u(j+1); % angular velocity command
        x_new(i) = x_old(i) + dt * v * cos(x_old(i+2));
        x_new(i+1) = x_old(i+1) + dt * v * sin(x_old(i+2));
        x_new(i+2) = x_old(i+2) + dt * w;
        j = j + 2;
    end
end

% generate Jacobian for dynamics
function A = DynamicsJacobian(x, u, dt)
    A = eye(length(x));
    j = 1;
    for i = 1:3:length(x)
        v = u(j);
        A(i,i+2) = -dt * v * sin(x(i+2));
        A(i+1,i+2) = dt * v * cos(x(i+2));
        j = j + 2;
    end
end

% nonlinear measurement (range and bearing)
function y = g(x)
    y = zeros(length(x), 1);
    for i = 1:3:length(x)
        rho = norm(x(i:i+1)); % range
        unit = x(i:i+1) / rho; % bearing
        y(i) = rho; 
        y(i+1:i+2) = unit;
    end
end

% generate Jacobian for measurements
function C = MeasurementJacobian(x)
    C = zeros(length(x));
    for i = 1:3:length(x)
        p = x(i:i+1); % extract position
        pos_norm = norm(p); % TODO: Guard against pos_norm=0 -> division by zero
        C_row1 = [p(1) / pos_norm, p(2) / pos_norm, 0]; % range
        C_rows23 = [-p * p' / pos_norm^3 + eye(2) / pos_norm, zeros(2,1)]; % bearing (unit vector)
        C(i,i:i+2) = C_row1;
        C(i+1:i+2,i:i+2) = C_rows23;
    end
end

