format longg
%% Task 1B
% Constants
tau_s = 1; % Propagation interval
S_a = 5; % PSD
Sigma = 2.5; % Measurement error standard deviation
% Extract data from table
Pos_ECEF_table = readtable("Workshop2_GNSS_Pos_ECEF.csv");
Pos_ECEF_arr = table2array(Pos_ECEF_table);
Pos_ECEF = Pos_ECEF_arr(:,2:end);
t = Pos_ECEF_arr(:,1);
%% Initial state estimation
x_est = [2447019;-5884199;-284783;184;77;0];
P_est = [100*eye(3), zeros(3);
            zeros(3), 25*eye(3)];
%% Transition matrix
Phi = [eye(3), tau_s*eye(3);
     zeros(3), eye(3)];
%% System noise covariance matrix
Q = [S_a*tau_s^3*eye(3)/3,  S_a*tau_s^2*eye(3)/2;
     S_a*tau_s^2*eye(3)/2,  S_a*tau_s*eye(3)];
%% Measurement matrix
H = [eye(3),zeros(3)];
%% Measurement noise covariance matrix
R = Sigma^2*eye(3);
%% Filter iteration
for i = 1:length(t)
    % Propagate state estimate
    x_p_est = Phi*x_est;
    % Propagate state error covariance matrix
    P_p_est = Phi*P_est*Phi' + Q;
    % Kalman gain matrix
    K = P_p_est*H'/(H*P_p_est*H'+R);
    % GNSS position solution
    r_ea = Pos_ECEF(i,:)';
    % Propagated position solution
    r_ea_pre = x_p_est(1:3,:);
    % Measurement innovation vector
    epsilon_z_pre = r_ea - r_ea_pre;
    % Update state estimation
    x_u_est = x_p_est + K*epsilon_z_pre;
    % Update error covariance matrix
    P_u_est = (eye(6) - K*H)*P_p_est;
    % Convert to NED formate
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_u_est(1:3),x_u_est(4:6));
    L_b = rad2deg(L_b);
    lambda_b = rad2deg(lambda_b);
    x_NED(i,1:3) = [L_b,lambda_b,h_b];
    x_NED(i,4:6) = v_eb_n';
    % Use updated results as previous results for the next iteration
    x_est = x_u_est;
    P_est = P_u_est;
end
disp(x_NED)