%% Task3
format longg
Pseudo_ranges_table = readtable("Workshop1_Pseudo_ranges.csv");
Pseudo_ranges_arr = table2array(Pseudo_ranges_table);
num_s = Pseudo_ranges_arr(1,2:end); % length = 8
time = Pseudo_ranges_arr(2:end,1); % length = 11
Pseudo_ranges = Pseudo_ranges_arr(2:end, 2:end);
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
c = 299792458; % Speed of light in m/s
Sigma_Rho = 5; % measurement error standard deviation
T = 6; % Outlier detection threshold
outlier = [];
residual = 0;
residual_pre = 0;
r_ea = get_r_ea(); % 3*1
% disp(r_ea)
r_es = get_r_es(); % 8*3*11
% disp(r_es)
r_as = zeros(length(time), length(num_s)); % 11*8
C_I_e = ones(3,3, length(num_s));
epsilon_Rho = 0; % predicted receiver clock offset
Rho_as = Pseudo_ranges(1,:)';
p_a = zeros(length(time),length(r_ea)); % 11*3

%% Iterate through time
for i = 1:length(time)
    r_es_i = r_es(:,:,i);
    num_s = Pseudo_ranges_arr(1,2:end); % Reset num_s to 8
    Rho_as = Pseudo_ranges(i,:)';
    isoutlier = 1;
    num_o = 0;
    while isoutlier == 1
        isoutlier = 0;
        u_as = zeros(size(r_es_i)); % num_s*3
        r_as_i = zeros(1,length(num_s) - num_o);
        % Iterate through satellite
        for j = 1: length(num_s) - num_o
            % Calculagte r_as_i
            r_as_i(j) = sqrt((r_es_i(j,:)' - r_ea)' * (r_es_i(j,:)' - r_ea));
            C_I_e(:,:,j) = [1,                      omega_ie*r_as_i(j)/c,   0;...
                            -omega_ie*r_as_i(j)/c,  1,                      0;...
                            0,                      0,                      1];
            diff = C_I_e(:,:,j)*r_es_i(j,:)' - r_ea;
            r_as_i(j) = sqrt(diff' * diff); % 1*num_s
            % Calculate u_as
            u_as(j,:) = diff/r_as_i(j);
        end
        
        %% Formulate x_pre, epsilon_z_pre, and H_G
        x_pre = [r_ea;epsilon_Rho];
        epsilon_z_pre = Rho_as - r_as_i' - epsilon_Rho;
        H_G = [-u_as,ones(length(num_s)-num_o,1)];
        % disp(H_G)
        
        %% Compute the residuals vector
        v = (H_G / (H_G' * H_G) * H_G' - eye(length(num_s)-num_o))*epsilon_z_pre; % 8*1
        % Compute the residuals covariance matrix
        C_v = (eye(length(num_s)-num_o) - H_G / (H_G' * H_G) * H_G')*Sigma_Rho^2; % 8*8
        for k = 1:length(v)
            residual =  abs(v(k)) - sqrt(C_v(k,k))*T;
            % Locate the time and satellite with largest residule
            if residual > residual_pre
                idx = k;
                residual_pre = residual;
                outlier = [time(i),num_s(idx)];
                isoutlier = 1;
                disp(['outlier occurred to satellite ' num2str(num_s(idx)) ' at time ' num2str(time(i))]);
            end
        end
        if isoutlier == 1
            num_o = num_o + 1;
            % Remove the r_es of the outlier satellite
            r_es_i(idx,:) = [];
            % Remove current pseudo range of the outlier satellite
            Rho_as(idx) = [];
        end
    end

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
disp('Outliers has been removed, no outlier exists now')
disp(p_a)