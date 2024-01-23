%% Task1 c) and d)
format longg
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
c = 299792458; % Speed of light in m/s
r_ea = zeros(3,1);
% disp(r_ea) % Task1 a)
r_es = get_r_es(); % 8*3*11
r_es_0 = r_es(:,:,1); % 8*3
% disp(r_es) % Task1 b)
r_as = zeros(length(r_es_0),1);
C_I_e = ones(3,3, length(r_es_0));
u_as = zeros(size(r_es_0));
epsilon_Rho = 0; % predicted receiver clock offset
Pseudo_ranges = readtable("Workshop1_Pseudo_ranges.csv");
Rho_as = table2array(Pseudo_ranges(2,2:end))';
x_est_a = [-4648608.2937218; 2555877.78608783; -3529219.8048316];
d = 100;
cnt = 0;

while d > 0.1
    for i = 1: length(r_es_0)
        r_as(i,:) = sqrt((r_es_0(i,:)' - r_ea)' * (r_es_0(i,:)' - r_ea));
        C_I_e(:,:,i) = [1,                      omega_ie*r_as(i,:)/c,   0;...
                        -omega_ie*r_as(i,:)/c,  1,                      0;...
                        0,                      0,                      1];
        diff = (C_I_e(:,:,i)*r_es_0(i,:)' - r_ea);
        r_as(i,:) = sqrt(diff' * diff);
        u_as(i,:) = (diff/r_as(i,:))';
    end
    % disp(r_as) % Task1 c)
    % disp(u_as) % Task1 d)
    
    %% Task 1 e)
    x_pre = [r_ea;epsilon_Rho];
    epsilon_z_pre = Rho_as - r_as - epsilon_Rho;
    H_G = [-u_as,ones(length(r_es_0),1)];
    % disp(H_G) % Task1 e)
    
    %% Task 1 f)
    x_est = x_pre + ((H_G'*H_G)\H_G')*epsilon_z_pre;
    % disp(x_est(1:3,:)) % Task1 f)
    x_est_b = x_est(1:3,:);

    %% Task 1 g)
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_est(1:3,:),0);
    L_b = rad2deg(L_b);
    lambda_b = rad2deg(lambda_b);
    p_a_b = [L_b;lambda_b;h_b];
    % disp(p_a_b) % Task1 g)
    d = sqrt(sum(x_est_a - x_est_b).^2);
    r_ea = x_est_b;
    cnt = cnt + 1;
    disp(cnt)
    disp(d)
end