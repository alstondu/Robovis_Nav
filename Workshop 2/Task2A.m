format longg
% Task 2A

%% Constants
tau_s = 1; % Propagation interval
S_a = 5; % acceleration power pectral PSD
S_rho = 10; % clock phase PSD
S_r = 0.05; % clock frequency PSD
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]); % Skew symmetric matrix of omega_ie
c = 299792458; % Speed of light in m/s

%% Initialization
[r_es, v_es] = get_x_es(); % 10*3*time
r_es_0 = r_es(:,:,1); % 10*3
v_es_0 = v_es(:,:,1); % 10*3
r_as = zeros(length(r_es_0),1); % 10*1
rate_as = zeros(size(r_as)); % 10*1
% Line of sight unit vector
u_as = zeros(size(r_es_0)); % 10*3
% Measurement innovation vector
epsilon_z_pre = zeros(20,1);
% Initial state estimation
[x_est,P_est] = Initialise_GNSS_KF;
% Sagnac effect compensation matrix
C_I_e = ones(3,3, length(r_es_0));  % 3*3*10

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

% Propagate state estimate
x_p_est = Phi*x_est;
% Propagate state error covariance matrix
P_p_est = Phi*P_est*Phi' + Q;

% Predicted Cartesian ECEF user position and velocity
r_ea = x_p_est(1:3,:);
v_ea = x_p_est(4:6,:);
% Calculate antenna-satellite range, line-of-signt unit vector, and predicted range rates
for i = 1: length(r_es_0)
    r_as(i,:) = sqrt((r_es_0(i,:)' - r_ea)' * (r_es_0(i,:)' - r_ea));
    C_I_e(:,:,i) = [1,                      omega_ie*r_as(i,:)/c,   0;...
                    -omega_ie*r_as(i,:)/c,  1,                      0;...
                    0,                      0,                      1];
    diff = (C_I_e(:,:,i)*r_es_0(i,:)' - r_ea);
    % Antenna-satellite range
    r_as(i,:) = sqrt(diff' * diff);
    % Line-of-signt unit vector
    u_as(i,:) = (diff/r_as(i,:))';
    % Predicted range rates
    rate_as(i,:) = u_as(i,:)*(C_I_e(:,:,i)*(v_es_0(i,:)' + Omega_ie*r_es_0(i,:)') - (v_ea + Omega_ie*r_ea));
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
% Measured pseudo-range
Pseudo_ranges = readtable("Workshop2_Pseudo_ranges.csv");
Rho_as = table2array(Pseudo_ranges(2,2:end))';
% Measured pseudo-range rate
Pseudo_rate = readtable("Workshop2_Pseudo_range_rates.csv");
Rho_as_dot = table2array(Pseudo_rate(2,2:end))';
% Measurement innovation vector
epsilon_z_pre(1:10) = Rho_as - r_as - epsilon_Rho;
epsilon_z_pre(11:20) = Rho_as_dot - rate_as - epsilon_Rho_dot;

% Update state estimation
x_u_est = x_p_est + K*epsilon_z_pre;
% Update error covariance matrix
P_u_est = (eye(8) - K*H)*P_p_est;

% Convert to NED formate
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_u_est(1:3),x_u_est(4:6));
L_b = rad2deg(L_b);
lambda_b = rad2deg(lambda_b);
x_NED = zeros(6,1);
x_NED(1:3) = [L_b;lambda_b;h_b];
x_NED(4:6) = v_eb_n;