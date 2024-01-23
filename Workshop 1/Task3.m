%% Task3
format longg
Pseudo_ranges_table = readtable("Workshop1_Pseudo_ranges.csv");
Pseudo_ranges = table2array(Pseudo_ranges_table);
num_s = Pseudo_ranges(1,2:end); % length = 8
time = Pseudo_ranges(2:end,1); % length = 11
Pseudo_ranges = Pseudo_ranges(2:end, 2:end);
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
c = 299792458; % Speed of light in m/s
r_ea = get_r_ea(); % 3*1
% disp(r_ea)
r_es = get_r_es(); % 8*3*11
% disp(r_es)
r_as = zeros(length(time), length(num_s)); % 11*8
C_I_e = ones(3,3, length(num_s));
u_as = zeros(size(r_es(:,:,1)));
epsilon_Rho = 0; % predicted receiver clock offset
Rho_as = Pseudo_ranges(1,:)';
p_a = zeros(length(time),length(r_ea)); % 11*3

%% Iterate through time
for i = 1:length(time)
    Rho_as = Pseudo_ranges(i,:)';
    % Iterate through satellite
    for j = 1: length(num_s) % 1-8
        % Calculagte r_as
        r_as(i,j) = sqrt((r_es(j,:,i)' - r_ea)' * (r_es(j,:,i)' - r_ea));
        C_I_e(:,:,j) = [1,                      omega_ie*r_as(i,j)/c,   0;...
                        -omega_ie*r_as(i,j)/c,  1,                      0;...
                        0,                      0,                      1];
        diff = C_I_e(:,:,j)*r_es(j,:,i)' - r_ea;
        r_as(i,j) = sqrt(diff' * diff);
        % Calculate u_as
        u_as(j,:) = diff/r_as(i,j);
    end
    % disp(r_as)
    % disp(u_as)
    
    %% Formulate x_pre, epsilon_z_pre, and H_G
    x_pre = [r_ea;epsilon_Rho];
    epsilon_z_pre = Rho_as - r_as(i,:)' - epsilon_Rho;
    H_G = [-u_as,ones(length(num_s),1)];
    % disp(H_G) % Task1 e)
    
    %% Calculate x_est
    x_est = x_pre + ((H_G'*H_G)\H_G')*epsilon_z_pre;
    % disp(x_est(1:3,:))
    r_ea_est = x_est(1:3,:);

    %% Convert x_est to latitude, longitude, and height
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_ea_est,0);
    L_b = rad2deg(L_b);
    lambda_b = rad2deg(lambda_b);
    p_a(i,:) = [L_b,lambda_b,h_b];

    %% Update location and receiver clock offset
    r_ea = r_ea_est;
    epsilon_Rho = x_est(4,:);
end

disp(p_a)