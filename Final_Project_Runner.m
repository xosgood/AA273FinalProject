%Final Project Main Runner

%Conditions that do not change for any simulation
v_follower_thresh = 5; %m/s for max speed difference in a step
omega_follower_thresh = 2.8648; %max degree difference in a step
const_data = struct();
const_data.n_F = 3; %Number Follower  States
const_data.K_p1 = [0.001; 0.001]; % Good Controller Gains
const_data.K_p2 = [0.1; 0.01]; % Bad Controller Gains 
const_data.K_d = [1;1]; %Derivative Controller Gains
const_data.t_delay = 0.3;
const_data.threshold = [v_follower_thresh; omega_follower_thresh];

%Conditions that can change by simulation
dt = 0.1;
run1 = struct(); 
run1.num_follow = 4; %number of followers
run1.add_noise = true;
run1.Qtrue = eye(const_data.n_F)*0.01;
run1.QtrueL =  0.01 * dt * eye(const_data.n_F);
run1.t_f = 10; %[s] runtime
run1.K_cycle = 100; %[s] switch between two controller mechanisms
run1.rho_error = [10, 50]; %Can be 10% shorter or 100% longer
run1.ang_error = 5; %Degrees Initial Offset can bet +/- the normal angle
run1.head_error = 30; %Degress initial heading angle is off (0 would be no error)
run1.rho_nom = 4; % meters
run1.ang_nom = 34; %Degrees
run1.dt = dt; % timestep
%Leader Dynamics
v = 10 * ones(1, length(0:dt:run1.t_f)); % velocity command
omega = sin((0:dt:run1.t_f)/3); % angular velocity command zeros(1, length(0:dt:run1.t_f));
u_L = [v; omega];
run1.u_L = u_L;

run1.control_from = "leader"; %"leader"; %"follower"
%If leader: all control is based off of desired dynamics
%If follower: control based off actual position of the bird's sub-leader

ekf.R = 2 * diag([0.1, 0.01, 0.01]);
ekf.Q = 1 * eye(const_data.n_F);
 

objs = computeICs(run1, const_data, true);
objs = runSimulation(objs, run1.num_follow, run1.u_L, run1.add_noise);
plotTrajs(objs, run1)

run1.control_from = "follower";
objs2 = computeICs(run1, const_data, false);
objs2 = runSimulation(objs2, run1.num_follow, run1.u_L, run1.add_noise);
plotTrajs(objs2, run1)




function objs = runSimulation(objs, num_follow, u_L, add_noise)
    prior_state = zeros(objs(1).n_F, num_follow+1);
    
    for i = 2:objs(1).numsteps

        for j = 1:(num_follow+1)
            if i-1-objs(1).iter_delay >= 1
                prior_state(:, j) = objs(j).x_F_act(:, i-1-objs(1).iter_delay);
                a =5;
            end
            
            if j == 1
                override_control = [1, u_L(1, i), u_L(2, i)];
            else
                override_control = [0, 0, 0];
            end
            
            objs(j).curr_ind = i;
            objs(j).desiredDynamics(override_control); 
            objs(j).actualDynamics(override_control, add_noise); 
            objs(j).generateControl(override_control, prior_state);
        
        end
        
    end

end

function objs = computeICs(run, const_data, plotme)  
    run.K_cycle = run.K_cycle*(1/run.dt); %convert into proper time unit
    iter_delay = round(const_data.t_delay / run.dt);
    const_data.threshold(2) = deg2rad(const_data.threshold(2));
    run.ang_nom = deg2rad(run.ang_nom);
        
    icD = zeros(run.num_follow+1, const_data.n_F);
    icA = zeros(run.num_follow+1, const_data.n_F);
    
    rho_range = [(100-run.rho_error(1))/100, (100+run.rho_error(2))/100]*run.rho_nom;
    ang_range = [run.ang_nom+deg2rad(run.ang_error), run.ang_nom-deg2rad(run.ang_error)];   
    head_range = [deg2rad(run.head_error), -deg2rad(run.head_error)];
    
    %generating random rho and angle values for actual starting condition
    rhos = (rho_range(2)-rho_range(1)).*rand(run.num_follow+1, 1) + rho_range(1);
    angs = (ang_range(1)-ang_range(2)).*rand(run.num_follow+1, 1) + ang_range(2);
    heads = (head_range(1)-head_range(2)).*rand(run.num_follow+1, 1) + head_range(2);
    
    %populating initial conditions
    for inc = 1:(run.num_follow+1)
        if inc == 1
            icD(inc, :) = [0, 0, 0];
            icA(inc, :) = [0, 0, 0]; %leader bird doesn't have an actual state just a true one
        elseif inc == 2
            icD(inc, :) = [icD(inc-1, 1)-run.rho_nom*cos(run.ang_nom), icD(inc-1, 2)+run.rho_nom*sin(run.ang_nom), 0];
            icA(inc, :) = [icD(inc-1, 1)-rhos(inc)*cos(angs(inc)), icD(inc-1, 2)+rhos(inc)*sin(angs(inc)), heads(inc)];
        elseif mod(inc, 2) == 1
            icD(inc, :) = [icD(inc-2, 1)-run.rho_nom*cos(run.ang_nom), icD(inc-2, 2)-run.rho_nom*sin(run.ang_nom), 0];
            icA(inc, :) = [icD(inc-2, 1)-rhos(inc)*cos(angs(inc)), icD(inc-2, 2)-rhos(inc)*sin(angs(inc)), heads(inc)];
        elseif mod(inc, 2) == 0
            icD(inc, :) = [icD(inc-2, 1)-run.rho_nom*cos(run.ang_nom), icD(inc-2, 2)+run.rho_nom*sin(run.ang_nom), 0];
            icA(inc, :) = [icD(inc-2, 1)-rhos(inc)*cos(angs(inc)), icD(inc-2, 2)+rhos(inc)*sin(angs(inc)), heads(inc)];
        else
            disp("there is an error somewhere");
        end
    end
    if(plotme)
        sz = 25;
        if mod(run.num_follow, 2) == 1
            odd_end = length(icD(:, 1))-1;
            even_end = length(icD(:, 1));
        else
            odd_end = length(icD(:, 1));
            even_end = length(icD(:, 1))-1;
        end
        figure;
        hold on; grid on; grid minor;
        scatter(icD(1, 1), icD(1,2), sz, 'g', 'filled')
        scatter(icD(2:end, 1), icD(2:end, 2), sz, 'b', 'filled')
        scatter(icA(2:end, 1), icA(2:end, 2), sz, 'r', 'filled');
        plot(icD((1:2:odd_end), 1), icD((1:2:odd_end), 2), 'k--', 'LineWidth', 1.5)
        plot(icA((1:2:odd_end), 1), icA([1:2:odd_end], 2), '--', 'Color', [105,105,105]/255, 'LineWidth', 1.5)
        plot(icD([1,2:2:even_end], 1), icD([1,2:2:even_end], 2), 'k--', 'LineWidth', 1.5)
        plot(icA([1,2:2:even_end], 1), icA([1,2:2:even_end], 2), '--', 'color', [105,105,105]/255, 'LineWidth', 1.5)
        scatter(icD(1, 1), icD(1,2), sz, 'g', 'filled')
        scatter(icD(2:end, 1), icD(2:end, 2), sz, 'b', 'filled')
        scatter(icA(2:end, 1), icA(2:end, 2), sz, 'r', 'filled');
        xlabel("x")
        ylabel("y")
        legend("Leader", "Follower Desired", "Follower Actual", "Desired Formation", "Actual Formation");
        hold off;
    end
    objs = Follower.empty;
    for i = 1:(run.num_follow+1)
        if i == 1
            Qtrue = run.QtrueL;
        else
            Qtrue = run.Qtrue;
        end
        objs(i) = Follower(const_data.n_F, run.t_f, run.dt, const_data.threshold, run.K_cycle, const_data.K_p1, const_data.K_p2, const_data.K_d, iter_delay, icD(i, :), icA(i, :), run.control_from, Qtrue, i);
    
    end
    
end

function createEKS(objs, ekf, run, const_data)
    
    ekf.n_S = run.num_follow * const_data.n_F;
    ekf.nl_dyns = zeros(ekf.n_S, objs(1).numsteps);
    ekf.mu_prior = zeros(ekf.n_S, objs(1).numsteps-1);
    ekf.mu_posterior = zeros(ekf.n_S, objs(1).numsteps-1);
    ekf.sigma_prior = zeros(ekf.n_S, ekf.n_S, objs(1).numsteps);
    ekf.sigma_posterior = zeros(ekf.n_S, ekf.n_S, objs(1).numsteps);
    if strcmpi(MEAS, "BOTH")
        ekf.nummeas = run.num_follow*3;
    elseif strcmpi(MEAS, "BEARING")
        ekf.nummeas = run.num_follow*2;
    elseif strcmpi(MEAS, "RANGE")
        ekf.nummeas = run.num_follow;
    else
        error("Error on Measurement Selection Type");
    end
    ekf.y = zeros(ekf.nummeas, objs(1).numsteps); % measurements
    ekf.y_hat = zeros(ekf.nummeas, objs(1).numsteps);
    
    ekf.mu_posterior(:,1) = -ones(ekf.n_S, 1);
    ekf.sigma_posterior(:,:,1) = eye(n_S);
    
    %Getting the N-L Dynamics Model
    indices = 1:1:const_data.n_F;
    for i = 1:run.num_follow
        ekf.nl_dyns(indices, :) = objs(i).x_F_act(:, :);
        indices = indices + const_data.n_F;
    end
    
    for i = 2:objs(1).numsteps
        x = ekf.nl_dyns(:, i);
        y_all_meas = g(x);
        %mu_prior(:, i-1)
        
        if strcmpi(MEAS, "BOTH")
            ekf.y(:, i) = y_all_meas;
        elseif strcmpi(MEAS, "BEARING")
            k = 2;
            states_pull = [];
            for j = 1:run.num_follow
                states_pull = [states_pull; k, k+1];
                k = k+3;
            end
            ekf.y(:, i) = y_all_meas(states_pull);
        elseif strcmpi(MEAS, "RANGE")
            states_pull = [];
            k = 1;
            for j = 1:run.num_follow
                states_pull = [states_pull; k];
                k = k+3;
            end
            ekf.y(:, i) = y_all_meas(states_pull);
        end
        
        

    end

end

function plotTrajs(objs, run)
    colors = getColors();
    [h, ~] = size(colors);
    if run.num_follow+1 > h
        num_rows = run.num_follow+1 - h;
        extra_colors = rand(num_rows, 3);
        colors = [colors; extra_colors];
    end
    x_L = objs(1).x_F_des(1:2, :);
    figure;
    hold on; grid on; grid minor;
    for i = 1:(run.num_follow+1)
        if i == 1
            str1 = "Leader";
            plot(x_L(1, :), x_L(2,:), 'b-', 'DisplayName', char(str1), 'LineWidth', 1.5, 'Color', colors(i, :))
        else
            str1 = "Desired Follower " + num2str(i-1);
            str2 = "Actual Follower  " + num2str(i-1);
            plot(x_L(1,:) + objs(i).x_F_des(1, :), x_L(2,:) + objs(i).x_F_des(2, :), ':', 'LineWidth', 1.5, 'DisplayName', char(str1), 'Color', colors(i, :))
            plot(x_L(1,:) + objs(i).x_F_act(1, :), x_L(2,:) + objs(i).x_F_act(2, :), 'o', 'MarkerSize', 5, 'DisplayName', char(str2), 'Color', colors(i, :))
        end
        drawnow;
    end
    
    hold off;



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

function C = MeasurementJacobian1(x)
    global MEAS;
    
    p = x(1:2); % extract position
    pos_norm = norm(p); % TODO: Guard against pos_norm=0 -> division by zero
    
    if MEAS == "RANGE"
        C_row1 = [p(1) / pos_norm, p(2) / pos_norm, 0]; % range
        C = C_row1;
    elseif MEAS == "BEARING"
        C_rows23 = [-p * p' / pos_norm^3 + eye(2) / pos_norm, zeros(2,1)]; % bearing (unit vector)
        C = C_rows23;
    elseif MEAS == "BOTH"
        C_row1 = [p(1) / pos_norm, p(2) / pos_norm, 0]; % range
        C_rows23 = [-p * p' / pos_norm^3 + eye(2) / pos_norm, zeros(2,1)]; % bearing (unit vector)
        C = [C_row1; C_rows23];
    else
        error("Not allowed 'MEAS' value")
    end
end

%returns a set of colors that are high contrast
function colors = getColors()
colors =    [[1 0 0];...
            [0 1 0];...
            [0 0 1];...
            [0 1 1];...
            [1 0 1];...
            [0 0 0];...
            [34 49 29]/255;...
            [0 0.4470 0.7410];...
            [0.8500 0.3250 0.0980];...
            [0.9290 0.6940 0.1250];...
            [0.4940 0.1840 0.5560];...
            [0.4660 0.6740 0.1880];...
            [0.3010 0.7450 0.9330];...
            [0.6350 0.0780 0.1840]];

end



