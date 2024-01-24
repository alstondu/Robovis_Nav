format longg
%% Task 1A
% Constants
tau_s = 1; % Propagation interval
S_a = 5; % PSD
Sigma = 2.5; % Measurement error standard deviation
Pos_ECEF_table = readtable("Workshop2_GNSS_Pos_ECEF.csv");
Pos_ECEF_arr = table2array(Pos_ECEF_table);
Pos_ECEF = Pos_ECEF_arr(:,2:end);
%% Task 1A a)
x0_est = [2447019;-5884199;-284783;184;77;0];
P0_est = [100*eye(3), zeros(3);
            zeros(3), 25*eye(3)];
%% Task 1A b)
Phi = [eye(3), tau_s*eye(3);
     zeros(3), eye(3)];
%% Task 1A c)
Q = [S_a*tau_s^3*eye(3)/3,  S_a*tau_s^2*eye(3)/2;
     S_a*tau_s^2*eye(3)/2,  S_a*tau_s*eye(3)];
%% Task 1A d)
x_p_est = Phi*x0_est;
%% Task 1A e)
P_p_est = Phi*P0_est*Phi' + Q;
%% Task 1A f)
H = [eye(3),zeros(3)];
%% Task 1A g)
R = Sigma^2*eye(3);
%% Task 1A h)
K = P_p_est*H'/(H*P_p_est*H'+R);
%% Task 1A i)
r_ea = Pos_ECEF(1,:)';
r_ea_pre = x_p_est(1:3,:);
epsilon_z_pre = r_ea - r_ea_pre;
%% Task 1A j)
x_u_est = x_p_est + K*epsilon_z_pre;
%% Task 1A k)
P_u_est = (eye(6) - K*H)*P_p_est;
%% Task 1A l)
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_u_est(1:3),x_u_est(4:6));
L_b = rad2deg(L_b);
lambda_b = rad2deg(lambda_b);
x_NED = zeros(6,1);
x_NED(1:3) = [L_b;lambda_b;h_b];
x_NED(4:6) = v_eb_n;