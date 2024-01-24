format longg
% Task 2B

%% Constants
tau_s = 1; % Propagation interval
S_a = 5; % acceleration power pectral PSD
S_rho = 10; % clock phase PSD
S_r = 0.05; % clock frequency PSD
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]); % Skew symmetric matrix of omega_ie
c = 299792458; % Speed of light in m/s

%% Extract data from tables
Pseudo_ranges_table = readtable("Workshop2_Pseudo_ranges.csv");
Pseudo_ranges = table2array(Pseudo_ranges_table);
Pseudo_rate_table = readtable("Workshop2_Pseudo_range_rates.csv");
Pseudo_rate = table2array(Pseudo_rate_table);
num_s = Pseudo_ranges(1,2:end); % length = 10
t = Pseudo_ranges(2:end,1); % length = 181
Pseudo_ranges = Pseudo_ranges(2:end, 2:end);
Pseudo_rate = Pseudo_rate(2:end, 2:end);

%% Initialization
[r_es, v_es] = get_x_es(); % 10*3*time
r_as = zeros(length(t), length(num_s)); % 181*10
rate_as = zeros(length(t), length(num_s)); % 181*10
% Line of sight unit vector
u_as = zeros(length(num_s),3); % 10*3
% Measurement innovation vector
epsilon_z_pre = zeros(20,1);
% Initial state estimation
[x_est,P_est] = Initialise_GNSS_KF;
% Sagnac effect compensation matrix
C_I_e = ones(3,3, length(num_s));  % 3*3*10
% Initialize measured pseudo-range
Rho_as = Pseudo_ranges(1,:)';
% Initialize measured pseudo-range rate
Rho_as_dot = Pseudo_rate(1,:)';
x_NED = zeros(length(t),6);

% Transition matrix
Phi = [eye(3), tau_s*eye(3), zeros(3,1), zeros(3,1);
     zeros(3),       eye(3), zeros(3,1), zeros(3,1);
   zeros(1,3),   zeros(1,3),          1,      tau_s;
   zeros(1,3),   zeros(1,3),          0,          1];

% System noise covariance matrix
Q = [S_a*tau_s^3*eye(3)/3,  S_a*tau_s^2*eye(3)/2, zeros(3,1), zeros(3,1);
     S_a*tau_s^2*eye(3)/2,      S_a*tau_s*eye(3), zeros(3,1), zeros(3,1);
               zeros(1,3),            zeros(1,3), S_cp*tau_s + S_cf*tau_s^3/3, S_cf*tau_s^2/2;
               zeros(1,3),            zeros(1,3), S_cf*tau_s^2/2,              S_cf*tau_s];
%% Iteration
% Iterate throught time
for i = 1:length(t)
    % Extract current data
    Rho_as = Pseudo_ranges(i,:)';
    Rho_as_dot = Pseudo_rate(i,:)';
    % Propagate state estimate
    x_p_est = Phi*x_est;
    % Propagate state error covariance matrix
    P_p_est = Phi*P_est*Phi' + Q;
    
    % Predicted Cartesian ECEF user position and velocity
    r_ea = x_p_est(1:3,:);
    v_ea = x_p_est(4:6,:);

    % Calculate antenna-satellite range, line-of-signt unit vector, and predicted range rates
    % Iterate through satellite
    for j = 1: length(num_s) % 1-10
        r_as(i,j) = sqrt((r_es(j,:,i)' - r_ea)' * (r_es(j,:,i)' - r_ea));
        C_I_e(:,:,j) = [1,                      omega_ie*r_as(i,j)/c,   0;...
                        -omega_ie*r_as(i,j)/c,  1,                      0;...
                        0,                      0,                      1];
        diff = C_I_e(:,:,j)*r_es(j,:,i)' - r_ea;
        % Antenna-satellite range
        r_as(i,j) = sqrt(diff' * diff);
        % Line-of-signt unit vector
        u_as(j,:) = diff/r_as(i,j);
        % Predicted range rates
        rate_as(i,j) = u_as(j,:)*(C_I_e(:,:,j)*(v_es(j,:,i)' + Omega_ie*r_es(j,:,i)') - (v_ea + Omega_ie*r_ea));
    end
    % Measurement matrix
    H = [-u_as,              zeros(size(u_as)),  ones(10,1),   zeros(10,1);
        zeros(size(u_as)),   -u_as,              zeros(10,1),  ones(10,1)];
    
    % Measurement noise covariance matrix
    R = [S_rho^2*eye(10,10), zeros(10,10);
               zeros(10,10), S_r^2*eye(10,10)];
    
    % Kalman gain matrix
    K = P_p_est*H'/(H*P_p_est*H'+R);
    
    % propagated receiver clock offset
    epsilon_Rho = x_p_est(7); 
    % propagated receiver clock drift
    epsilon_Rho_dot = x_p_est(8);
    
    % Measurement innovation vector
    epsilon_z_pre(1:10) = Rho_as - r_as(i,:)' - epsilon_Rho;
    epsilon_z_pre(11:20) = Rho_as_dot - rate_as(i,:)' - epsilon_Rho_dot;
    
    % Update state estimation
    x_u_est = x_p_est + K*epsilon_z_pre;
    % Update error covariance matrix
    P_u_est = (eye(8) - K*H)*P_p_est;
    
    % Convert to NED formate
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_u_est(1:3),x_u_est(4:6));
    L_b = rad2deg(L_b);
    lambda_b = rad2deg(lambda_b);
    x_NED(i,1:3) = [L_b,lambda_b,h_b];
    x_NED(i,4:6) = v_eb_n';

    % Update states and covariance for iteration
    x_est = x_u_est;
    P_est = P_u_est;
end