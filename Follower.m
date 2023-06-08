%Follower Dynamics Class
classdef Follower < handle

    properties 
        init_condD  double {mustBeNonNan} = [0; 0; 0;] % Desired x, y, heading angle
        init_condA  double {mustBeNonNan} = [0; 0; 0;] % Desired x, y, heading angle
        dt          double {mustBePositive} = 0.1      % time step
        x_F_act     % Actual Follower State Values (Noisy and Full of Error)
        x_F_des     % Desired Follower State Values (Essentially Offset and Pinned to Leader Dynamics
        u_F         % Controller Commands
        tspan       % Time Span
        n_F         double {mustBeNonnegative} = 3 %Number of States
        numsteps    % Number of Time Steps
        curr_ind    double {mustBeNonnegative} = 1
        Qtrue       % True Noise Dynamics
        K_cycle     %how long to stay in good or bad controller gains [s]
        K_p1        %good controller gains
        K_p2        %bad controller gains
        K_d         % Derivative Controller Gains
        v_follower_max_thresh double {mustBeNonnegative} = 5; % Max Velocity Change Per Time Step
        omega_follower_max_thresh double {mustBeNonnegative} = 0.05; %Max Angle Change Per Time Step
        e_prev      %Previous Error Signal
        iter_delay  double {mustBeNonnegative} = 3
        controlFrom string {mustBeText}
        birdnum    % bird number
    end
    
    methods
        function obj = Follower(num_states, tf, dt, threshold, k_cycle, k_p1, k_p2, k_d, iter_delay, initial_condD, initial_condA, control_from, Qtrue, bird_num)
            %Constructor for Follower Class
            if nargin < 14
                disp(nargin)
                error("Incorrect Number of Input Entries.")
            else
                obj.dt = dt;
                obj.tspan = 0:dt:tf;
                obj.numsteps = length(obj.tspan);
                obj.n_F = num_states;
                obj.x_F_act = zeros(obj.n_F, obj.numsteps);
                obj.x_F_des = zeros(obj.n_F, obj.numsteps);
                obj.v_follower_max_thresh = threshold(1);
                obj.omega_follower_max_thresh = threshold(2);
                obj.u_F = zeros(2, obj.numsteps);
                obj.K_cycle = k_cycle;
                obj.K_p1 = k_p1;
                obj.K_p2 = k_p2;
                obj.K_d = k_d;
                obj.init_condD = initial_condD;
                obj.init_condA = initial_condA;
                obj.controlFrom = control_from;
                obj.Qtrue = Qtrue;
                obj.iter_delay = iter_delay;
                obj.birdnum = bird_num;
            end
            obj.x_F_des(:, 1) = obj.init_condD;
            obj.x_F_act(:, 1) = obj.init_condA;
            obj.e_prev = zeros(obj.n_F,1);
        end
         
        function desiredDynamics(obj, override_control)   
            t = obj.curr_ind;
            if boolean(override_control(1)) %for leader control
                u = [override_control(2); override_control(3)];
            else
                u = [0; 0]; 
            end
            obj.x_F_des(1, t) = obj.x_F_des(1, t-1) + obj.dt * u(1) * cos(obj.x_F_des(3, t-1));
            obj.x_F_des(2, t) = obj.x_F_des(2, t-1) + obj.dt * u(1) * sin(obj.x_F_des(3, t-1));
            obj.x_F_des(3, t) = obj.x_F_des(3, t-1) + obj.dt * u(2);
        end
        
        function actualDynamics(obj, override_control, add_noise)
            t = obj.curr_ind;
            if boolean(override_control(1))  %for leader control
                u = [override_control(2); override_control(3)];
                w_rel = zeros(obj.n_F, 1); %leader does not have noise
            else
                u = obj.u_F(:, t-1);
                if add_noise
                    w_rel = mvnrnd(zeros(obj.n_F,1), obj.Qtrue)';
                else
                    w_rel = mvnrnd(zeros(obj.n_F,1));
                end
            end
            obj.x_F_act(1, t) = obj.x_F_act(1, t-1) + obj.dt * u(1) * cos(obj.x_F_act(3, t-1));
            obj.x_F_act(2, t) = obj.x_F_act(2, t-1) + obj.dt * u(1) * sin(obj.x_F_act(3, t-1));
            obj.x_F_act(3, t) = obj.x_F_act(3, t-1) + obj.dt * u(2);
            obj.x_F_act(:, t) = obj.x_F_act(:, t) + w_rel;
        end
            
        function generateControl(obj, override_control, prior_state)
            disp(" ");
            i = obj.curr_ind;
            if boolean(override_control(1)) %for leader control or something else
                obj.u_F(:, i) = [override_control(2); override_control(3)];
                %disp(['Bird Num: ', num2str(obj.birdnum), ' Control Overridden']);
                %disp(['Control Output u_F1: ', num2str(obj.u_F(1,i)), ' u_F2: ', num2str(obj.u_F(2,i))]);
            else
                if strcmpi(obj.controlFrom, "leader")
                    %disp(['Bird Num: ', num2str(obj.birdnum), ' Inside Track Leader Bird']);
                    % control law for follower (control for next timestep)
                    if i-1-obj.iter_delay >= 1
                        x_F_actT = obj.x_F_des(:,i-1-obj.iter_delay);
                    else
                        x_F_actT = zeros(obj.n_F,1);
                    end
                    obj.u_F(:, i) = implementControl(obj, i, x_F_actT, "Desired");
                else
                    %disp(['Bird Num: ', num2str(obj.birdnum), ' Inside Track Sub Leader Bird']);
                    if obj.birdnum == 1 || obj.birdnum == 2 || obj.birdnum == 3
                        if i-1-obj.iter_delay >= 1
                            x_F_actT = obj.x_F_des(:,i-1-obj.iter_delay);
                        else
                            x_F_actT = zeros(obj.n_F,1);
                        end
                        obj.u_F(:, i) = implementControl(obj, i, x_F_actT, "Desired");
                    else
                        x_F_actT = prior_state(:, obj.birdnum-2);
                        obj.u_F(:, i) = implementControl(obj, i, x_F_actT, "Actual");
                        if i-1-obj.iter_delay >= 1
                            x_F_test = obj.x_F_des(:,i-1-obj.iter_delay);
                        else
                            x_F_test = zeros(obj.n_F,1);
                        end
                        
                        aaa = implementControl(obj, i, x_F_test, "Specird");
                    end
                end
            end
        end
        
        function u_F = implementControl(obj, i, x_F_actT, str1)
            u_F = zeros(length(obj.u_F(:,1)),1);
            if i-1-obj.iter_delay >= 1
                e = obj.x_F_act(:,i-1-obj.iter_delay) - x_F_actT(:);
                %disp(['Error Used e1: ', num2str(e(1)), ' e2: ', num2str(e(2)), ' e3: ', num2str(e(3))]);
                e11 = obj.x_F_act(:,i-1-obj.iter_delay) - obj.x_F_des(:,i-1-obj.iter_delay);
                %disp(['Error Lead e1: ', num2str(e11(1)), ' e2: ', num2str(e11(2)), ' e3: ', num2str(e11(3))]);
            else
                e = obj.e_prev;
            end
            % P control
            e_x = e(1);
            e_y = e(2);
            e_psi = e(3);
            e_rho = norm(e(1:2));
            Beta = atan2(e_y, e_x);
            if mod(i, 2 * obj.K_cycle) <= obj.K_cycle
                K_p = obj.K_p1;
            else
                K_p = obj.K_p2;
            end

            % D control
            e_diff = e - obj.e_prev;
            obj.e_prev = e; % store e_prev for derivative control
            e_diff_x = e_diff(1);
            e_diff_y = e_diff(2);
            e_diff_psi = e_diff(3);
            e_diff_rho = norm(e_diff(1:2));
            %disp([char(str1), ' P e_Norm: ', num2str(e_rho), ' D e_Norm: ', num2str(e_diff_rho)]);
            Beta_diff = atan2(e_diff_y, e_diff_x);

            % combine PD control
            u_F(:) = K_p .* [e_rho;
                               e_psi - Beta] ...
                     + obj.K_d .* [e_diff_rho;
                               e_diff_psi - Beta_diff];

            u_F(1) = sign(u_F(1)) * min(abs(u_F(1)), obj.v_follower_max_thresh);
            u_F(2) = sign(u_F(2)) * min(abs(u_F(2)), obj.omega_follower_max_thresh);
            %disp([char(str1), ' Control Output u_F1: ', num2str(u_F(1)), ' u_F2: ', num2str(u_F(2))]);
        end
    end
end