format longg
% Task4

%% Extract data from tables
Pseudo_ranges_table = readtable("Workshop1_Pseudo_ranges.csv");
Pseudo_ranges = table2array(Pseudo_ranges_table);
Pseudo_rate_table = readtable("Workshop1_Pseudo_range_rates.csv");
Pseudo_rate = table2array(Pseudo_rate_table);
num_s = Pseudo_ranges(1,2:end); % length = 8
time = Pseudo_ranges(2:end,1); % length = 11
Pseudo_ranges = Pseudo_ranges(2:end, 2:end);
Pseudo_rate = Pseudo_rate(2:end, 2:end);

%% Constants
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]); % Skew symmetric matrix of omega_ie
c = 299792458; % Speed of light in m/s
epsilon_Rho = 0; % predicted receiver clock offset
epsilon_Rho_dot = 0; % predicted receiver clock drift

%% Initialization
r_ea = get_r_ea(); % 3*1
v_ea = zeros(size(r_ea));
% Initial state estimation
[r_es, v_es] = get_x_es(); % 8*3*11
r_as = zeros(length(time), length(num_s)); % 11*8
rate_as = zeros(size(r_as));
C_I_e = ones(3,3, length(num_s));
% Line of sight unit vector
u_as = zeros(size(r_es(:,:,1)));
% Initialize measured pseudo-range
Rho_as = Pseudo_ranges(1,:)';
% Initialize measured pseudo-range rate
Rho_as_dot = Pseudo_rate(1,:)';
v_a = zeros(length(time),length(r_ea)); % 11*3
p_a = zeros(length(time),length(r_ea)); % 11*3
% Measurement innovation vector
epsilon_z_pre = zeros(16,1);

%% Iteration
% Iterate through time
for i = 1:length(time)
    % Extract current data
    Rho_as = Pseudo_ranges(i,:)';
    Rho_as_dot = Pseudo_rate(i,:)';
    % Iterate through satellite
    for j = 1: length(num_s) % 1-8
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

    % Formulate x_pre, epsilon_z_pre, and H_G
    x_pre = [r_ea;v_ea;epsilon_Rho;epsilon_Rho_dot];
    epsilon_z_pre(1:8) = Rho_as - r_as(i,:)' - epsilon_Rho;
    epsilon_z_pre(9:16)= Rho_as_dot - rate_as(i,:)' - epsilon_Rho_dot;
    H_G = [-u_as,              zeros(size(u_as)),  ones(8,1),   zeros(8,1);
    zeros(size(u_as)),   -u_as,              zeros(8,1),  ones(8,1)];
    
    % Calculate x_est
    x_est = x_pre + ((H_G'*H_G)\H_G')*epsilon_z_pre;
    r_ea_est = x_est(1:3,:);
    v_ea_est = x_est(4:6,:);

    % Convertion from ECEF to NED
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_ea_est,v_ea_est);
    L_b = rad2deg(L_b);
    lambda_b = rad2deg(lambda_b);
    p_a(i,:) = [L_b,lambda_b,h_b];
    v_a(i,:) = v_eb_n;

    % Update velocity and receiver clock drift for iteration
    r_ea = r_ea_est;
    v_ea = v_ea_est;
    epsilon_Rho = x_est(7,:);
    epsilon_Rho_dot = x_est(8,:);
end