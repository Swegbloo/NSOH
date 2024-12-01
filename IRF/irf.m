% MATLAB code to calculate damping coefficients B(τ)
clc;
clear;
% Load the input matrix from a CSV file or other source
inputMatrix = load('data_am.csv'); 
forceMatrix = load('F_and_M.csv');
%% Part 1: MATLAB code to calculate the infinte added mass A(inf)
% axx = Axx, when frequency is infinite (for the highest frequency) 
A_33 = inputMatrix(46,6); 
A_35 = inputMatrix(46,7);
A_53 = inputMatrix(46,8);
A_55 = inputMatrix(46,9);

% Display the results
fprintf('Infinite Added Mass Coefficients:\n');
fprintf('A_33 = %.4f\n', A_33);
fprintf('A_35 = %.4f\n', A_35);
fprintf('A_53 = %.4f\n', A_53);
fprintf('A_55 = %.4f\n', A_55);

%% Part 2: Solve for damping term
% Load the input matrix from a CSV file or other source
% Example: Assume 'data.csv' contains the matrix with columns:
% Frequency | b_33 | b_35 | b_53 | b_55
inputMatrix = load('data.csv'); 

% Extract columns from the input matrix
frequency = inputMatrix(:, 1); % First column: frequency (ω)
b_33 = inputMatrix(:, 2);      % Second column: b_33 values
b_35 = inputMatrix(:, 3);      % Third column: b_35 values
b_53 = inputMatrix(:, 4);      % Fourth column: b_53 values
b_55 = inputMatrix(:, 5);      % Fifth column: b_55 values

% Load the τ matrix from a CSV file (assuming it is in a single column)
tau_values = readmatrix('tau_values.csv');  % Adjust file name if necessary

% Initialize variables
N = length(frequency); % Number of frequency rows
num_tau = length(tau_values); % Number of τ values

% Preallocate arrays for storing damping coefficients
B_33 = zeros(num_tau, 1);
B_35 = zeros(num_tau, 1);
B_53 = zeros(num_tau, 1);
B_55 = zeros(num_tau, 1);

% Loop through each τ value
for t = 1:num_tau
    tau = tau_values(t); % Current τ value
    
    % Initialize integrals for this τ
    I_33 = zeros(N, 1);
    I_35 = zeros(N, 1);
    I_53 = zeros(N, 1);
    I_55 = zeros(N, 1);
    
    % Compute integrals for each frequency interval
    for i = 2:N
        omega_prev = frequency(i-1); % Lower bound of frequency interval
        omega_curr = frequency(i);   % Upper bound of frequency interval
        
        % Extract b(ω) values for the current interval
        b33_vals = b_33(i-1:i);
        b35_vals = b_35(i-1:i);
        b53_vals = b_53(i-1:i);
        b55_vals = b_55(i-1:i);
        
        % Calculate the integral for each mode using interpolation
        I_33(i) = integral(@(omega) interp1([omega_prev, omega_curr], b33_vals, omega) .* cos(omega * tau), omega_prev, omega_curr);
        I_35(i) = integral(@(omega) interp1([omega_prev, omega_curr], b35_vals, omega) .* cos(omega * tau), omega_prev, omega_curr);
        I_53(i) = integral(@(omega) interp1([omega_prev, omega_curr], b53_vals, omega) .* cos(omega * tau), omega_prev, omega_curr);
        I_55(i) = integral(@(omega) interp1([omega_prev, omega_curr], b55_vals, omega) .* cos(omega * tau), omega_prev, omega_curr);
    end
    
    % Compute the damping coefficients for the current τ
    B_33(t) = (2 / pi) * sum(abs(I_33));
    B_35(t) = (2 / pi) * sum(abs(I_35));
    B_53(t) = (2 / pi) * sum(abs(I_53));
    B_55(t) = (2 / pi) * sum(abs(I_55));
end

% Combine τ and damping coefficients into a single matrix
result_matrix = [tau_values, B_33, B_35, B_53, B_55];

% Export results to a CSV file
writematrix(result_matrix, 'damping_coefficients.csv'); % Export to 'damping_coefficients.csv'

% Display results
fprintf('Damping Coefficients for each τ:\n');
for t = 1:num_tau
    fprintf('τ = %.2f: B_33 = %.4f, B_35 = %.4f, B_53 = %.4f, B_55 = %.4f\n', ...
        tau_values(t), B_33(t), B_35(t), B_53(t), B_55(t));
end
%% Part 3: solving for the motion

%initialize
x = zeros(num_tau,1);
v = zeros(num_tau,1);
a = zeros(num_tau,1);
xr = zeros(num_tau,1);
vr = zeros(num_tau,1);
ar = zeros(num_tau,1);
x(1,1) = 0;
v(1,1) = 0;
a(1,1) = 0;
xr(1,1) = 0;
vr(1,1) = 0;
ar(1,1) = 0;
damp_33 = 0;
damp_35 = 0;
damp_53 = 0;
damp_55 = 0;
dt = 0.3867;
F3 = forceMatrix(:,2);
F5 = forceMatrix(:,3);
M_33 = 22843.0 *1000;
Ky = 42.8;
I = M_33*Ky^2;
C_33 = 30448 *1000;
C_35 = -102927.337 *1000;
C_53 = -102927.337 *1000;
C_55 = 44366332 *1000;
%time marching algorithm
for i = 2:298
    for j = 1:297
        if i>j
            damp_33 = damp_33 + B_33(j)*v(i-j,1)*dt;
            damp_35 = damp_35 + B_35(j)*vr(i-j,1)*dt;
            damp_53 = damp_53 + B_53(j)*v(i-j,1)*dt;
            damp_55 = damp_55 + B_55(j)*vr(i-j,1)*dt;
        else
            break;
        end
    end
        M = [M_33+A_33 I+A_35;M_33+A_53 I+A_55];
        F = [F3(i)-C_33*x(i)-damp_33-C_35*xr(i)-damp_35;F5(i)-C_55*xr(i)-damp_55-C_53*x(i)-damp_53];
        X = M\F;
        a(i) = X(1,1);
        ar(i) = X(2,1);
        v(i) = v(i-1) + a(i-1)*dt;
        vr(i) = vr(i-1) + ar(i-1)*dt;
        x(i) = x(i-1) + v(i-1)*dt;
        xr(i) = xr(i-1) + vr(i-1)*dt;
end
hold on;
plot(tau_values,x);
plot(tau_values,xr);
hold off;