format longg
% Task 2

%% Constants
tau_s = 0.5; % Time interval
Sigma_v = 0.1; % DR velocity uncertainty
Sigma_r = 10; % DR position uncertainty
S_DR = 0.2; % DR velocity error PSD
Sigma_Gv = 0.02; % GNSS velocity uncertainty
Sigma_Gr = 5; % GNSS position uncertainty

%% Extract data
GNSS_data_table = readtable("Workshop3_GNSS_Pos_Vel_NED.csv");
GNSS_data = table2array(GNSS_data_table);
GNSS_data(:,2:3) = deg2rad(GNSS_data(:,2:3)); % Convert to rad
t = GNSS_data(:,1); % length = 351
% Load DR solution from Task1
load('DR_Solution.mat','DR_Solution');
DR_Solution(:,1:2) = deg2rad(DR_Solution(:,1:2)); % Convert to rad

%% Initialization
x = zeros(4,1); % States
L_k = GNSS_data(:,2); % Latitude
Lamda_k = GNSS_data(:,3); % Longitude
h = GNSS_data(:,4); % Height
v_NG = GNSS_data(:,5); %  GNSS velocity towards N
v_EG = GNSS_data(:,6); %  GNSS velocity towards E
GNSS_Solution = [L_k, Lamda_k, v_NG, v_EG];

% Initialize state estimation error covariance matrix
[R_N,R_E] = Radii_of_curvature(L_k(1));
P = [Sigma_v^2*eye(2), zeros(2);
    zeros(2),[Sigma_r^2/(R_N+h(1))^2, 0;
              0, Sigma_r^2/((R_E+h(1))^2*cos(L_k(1))^2)]];
% Initialize measurement matrix
H = [0,0,-1,0;
     0,0,0,-1;
    -1,0,0,0;
     0,-1,0,0];
epsilon_z_pre = zeros(4,1);
DR_Solution_C = zeros(size(DR_Solution));
DR_Solution_C(1,:) = DR_Solution(1,:); %？
% DR_Solution_C(1,:) = GNSS_Solution(1,:);

%% Iteration
for i = 2:length(t)
    % Current radii of curvatures
    % [R_N,R_E] = Radii_of_curvature(L_k(i-1));
    % Transition matrix
    Phi = [eye(2),zeros(2);
        [tau_s/(R_N + h(i-1)),0;
        0,tau_s/((R_E + h(i-1))*cos(L_k(i-1)))], eye(2)];
    % System noise covariance matrix
    Q = [S_DR*tau_s, 0, 0.5*S_DR*tau_s^2/(R_N + h(i-1)), 0;
        0, S_DR*tau_s, 0, 0.5*S_DR*tau_s^2/((R_E + h(i-1))*cos(L_k(i-1)));
        0.5*S_DR*tau_s^2/(R_N + h(i-1)), 0, (S_DR*tau_s^3/(R_N + h(i-1))^2)/3, 0;
        0, 0.5*S_DR*tau_s^2/((R_E + h(i-1))*cos(L_k(i-1))), 0, (S_DR*tau_s^3/((R_E + h(i-1))^2*cos(L_k(i-1))^2))/3];
    % Propagate state estimates
    x_p = Phi*x;
    % Propagate error covariance matrix
    P_p = Phi*P*Phi' + Q;
    % Update R_N and R_E
    [R_N,R_E] = Radii_of_curvature(L_k(i)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Measurement noise covariance matrix
    R = [[Sigma_Gr^2/(R_N + h(i)),0;
        0,Sigma_Gr^2/((R_E + h(i))^2*cos(L_k(i))^2)],zeros(2);
        zeros(2), Sigma_Gv^2*eye(2)];
    % Kalman gain matrix
    K = P_p*H'/(H*P_p*H'+R);
    % Measurement innovation vector
    epsilon_z_pre = GNSS_Solution(i,:)' - DR_Solution(i,:)' - H*x_p;
    % Update the state estimates
    x_u = x_p + K*epsilon_z_pre;
    % Update the error covariance matrix
    P_u = (eye(4) - K*H)*P_p;
    % Correct the DR solution with the Kalman filter
    DR_Solution_C(i,:) = DR_Solution(i,:)' - [x_u(3);x_u(4);x_u(1);x_u(2)];
    % Update states and covariance for iteration
    x = x_u;
    P = P_u;
end

% Convert the corrected solution back to degree
DR_Solution_C(:,1:2) = rad2deg(DR_Solution_C(:,1:2));
DR_Solution_C = [t, DR_Solution_C];
disp(DR_Solution_C)