format longg

% Task 1
%% Constants
L_k0 = deg2rad(50.4249580); % Initial Latitude
Lamda_k0 = deg2rad(-3.5957974); % Initial Longitude
h = 37.4; % Height
delta_t = 0.5; % Time interval

%% Extract data from tables
Speed_Heading_table = readtable("Workshop3_Speed_Heading.csv");
Speed_Heading = table2array(Speed_Heading_table);
t = Speed_Heading(:,1); % length = 351
v_m = Speed_Heading(:, 2); % Measured speed
Psi_m = Speed_Heading(:, 3); % Measured Heading in degree
Psi_m = deg2rad(Psi_m); % Measured Heading in radius

%% Initialization
v_a = zeros(2,1); % Average speed
L_k = zeros(length(t),1); % Latitude
L_k(1) = L_k0;
Lamda_k = zeros(length(t),1); % Longitude
Lamda_k(1) = Lamda_k0;
v_i = zeros(length(t),2); % Instantaneous velocity
v_i(1,1) = v_m(1)* cos(Psi_m(1));
v_i(1,2) = v_m(1)* sin(Psi_m(1));

%% Iteration
for i = 2:length(t)
    V_m2a = [cos(Psi_m(i)) + cos(Psi_m(i-1));
                sin(Psi_m(i)) + sin(Psi_m(i-1))];
    % Average speed
    v_a = 0.5*V_m2a*v_m(i);
    [R_N,R_E] = Radii_of_curvature(L_k(i-1));
    % latitude
    L_k(i) = L_k(i-1) + v_a(1)*delta_t/(R_N + h);
    % longitude
    Lamda_k(i) = Lamda_k(i-1) + v_a(2)*delta_t/((R_E + h)*cos(L_k(i)));
    % Instantaneous speed
    v_i(i,1) = 1.7*v_a(1) - 0.7*v_i(i-1,1);
    v_i(i,2) = 1.7*v_a(2) - 0.7*v_i(i-1,2);
end
L_k = rad2deg(L_k);
Lamda_k = rad2deg(Lamda_k);
DR_Solution = [L_k, Lamda_k, v_i];
% Save the data for Task 2
save('DR_Solution.mat','DR_Solution');
DR_Solution = [t,DR_Solution];
disp(DR_Solution)

