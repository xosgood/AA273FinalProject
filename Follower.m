%Follower Dynamics Class
classdef Follower < handle

    properties 
        init_condD double {mustBeNonNan} = [0; 0; 0;] % Desired x, y, heading angle
        init_condA double {mustBeNonNan} = [0; 0; 0;] % Desired x, y, heading angle
        dt         double {mustBePositive} = 0.1      % time step
        x_F_act    % Actual Follower State Values (Noisy and Full of Error)
        x_F_des    % Desired Follower State Values (Essentially Offset and Pinned to Leader Dynamics
        u_F        % Controller Commands
        tspan      % Time Span
        n_F        double {mustBeNonnegative} = 3 %Number of States
        numsteps   % Number of Time Steps
        curr_ind   double {mustBeNonnegative} = 1
        Qtrue      % True Noise Dynamics
        v_follower_max_thresh double {mustBeNonnegative} = 5; % Max Velocity Change Per Time Step
        omega_follower_max_thresh double {mustBeNonnegative} = 0.05; %Max Angle Change Per Time Step
        K_p 
    end
    
    methods
        function obj = Follower(num_states, tf, dt, initial_condD, initial_condA, Qtrue)
            %Constructor for Follower Class
            if nargin < 3
                error("Need to enter required number of follower state dims, final simulation time, and a time step.")
            else
                obj.dt = dt;
                obj.tspan = 0:dt:tf;
                obj.numsteps = length(obj.tspan);
                obj.n_F = num_states;
                obj.x_F_act = zeros(obj.n_F, obj.numsteps);
                obj.x_F_des = zeros(obj.n_F, obj.numsteps);
                
                obj.u_F = zeros(2, obj.numsteps);
                obj.K_p = [.1; .01] * 10 * sin(5*obj.tspan); % proportional gain for control loop for follower
            end
            
            if nargin == 3
                obj.init_condD = zeros(num_states, 1);
                obj.init_condA = zeros(num_states, 1);
                obj.Qtrue = eye(obj.n_F)*0.1;
            elseif nargin == 4
                obj.init_condD = initial_condD;
                obj.init_condA = zeros(num_states, 1);
                obj.Qtrue = eye(obj.n_F)*0.1;
            elseif nargin == 5
                obj.init_condD = initial_condD;
                obj.init_condA = initial_condA;
                obj.Qtrue = eye(obj.n_F)*0.1;
            else
                obj.init_condD = initial_condD;
                obj.init_condA = initial_condA;
                obj.Qtrue = Qtrue;
            end
            obj.x_F_des(:, 1) = obj.init_condD;
            obj.x_F_act(:, 1) = obj.init_condA;
        end
         
        function desiredDynamics(obj)   
            t = obj.curr_ind;
            u = [0; 0]; %obj.u_F(:, t);
            obj.x_F_des(1, t) = obj.x_F_des(1, t-1) + obj.dt * u(1) * cos(obj.x_F_des(3, t-1));
            obj.x_F_des(2, t) = obj.x_F_des(2, t-1) + obj.dt * u(1) * sin(obj.x_F_des(3, t-1));
            obj.x_F_des(3, t) = obj.x_F_des(3, t-1) + obj.dt * u(2);
        end
        
        function actualDynamics(obj, w_rel)
            t = obj.curr_ind;
            u = obj.u_F(:, t-1);
            obj.x_F_act(1, t) = obj.x_F_act(1, t-1) + obj.dt * u(1) * cos(obj.x_F_act(3, t-1));
            obj.x_F_act(2, t) = obj.x_F_act(2, t-1) + obj.dt * u(1) * sin(obj.x_F_act(3, t-1));
            obj.x_F_act(3, t) = obj.x_F_act(3, t-1) + obj.dt * u(2);
            obj.x_F_act(:, t) = obj.x_F_act(:, t) + w_rel;
        end
            
        function generateControl(obj)
            i = obj.curr_ind;
            u_F_temp = zeros(length(obj.x_F_act(:, i-1))-1,1);
            e = obj.x_F_act(:,i-1) - obj.x_F_des(:,i-1);
            e_x = e(1);
            e_y = e(2);
            e_psi = e(3);
            e_rho = norm(e(1:2));
            Beta = atan2(e_y, e_x);
            
            u_F_temp = obj.K_p(:,i) .* [e_rho;
                                    e_psi - Beta];
            u_F_temp(1) = sign(u_F_temp(1)) * min(abs(u_F_temp(1)), obj.v_follower_max_thresh);
            u_F_temp(2) = sign(u_F_temp(2)) * min(abs(u_F_temp(2)), obj.omega_follower_max_thresh);
            obj.u_F(:, i) = u_F_temp(:);
%             t = obj.curr_ind;
%             e = obj.x_F_act(:,t-1) - obj.x_F_des(:,t-1);
%             e_x = e(1);
%             e_y = e(2);
%             e_psi = e(3);
%             e_rho = norm(e(1:2));
%             Beta = atan2(e_y, e_x);
%             obj.u_F(:,t) = obj.K_p(:,t) .* [e_rho;
%                                     e_psi - Beta];
%             obj.u_F(1, t) = sign(obj.u_F(1, t)) * min(abs(obj.u_F(1, t)), obj.v_follower_max_thresh);
%             obj.u_F(2, t) = sign(obj.u_F(2, t)) * min(abs(obj.u_F(2, t)), obj.omega_follower_max_thresh);
        end
    end
end