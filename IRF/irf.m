clc;
clear;

% Load the input matrix from a CSV file or other source
inputMatrix = load('data.csv'); 

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


%% Part 2: MATLAB code to calculate damping coefficients B(τ)

% Example: Assume 'data.csv' contains the matrix with columns:
% Frequency | b_33 | b_35 | b_53 | b_55


% Extract columns from the input matrix
frequency = inputMatrix(:, 1); % First column: frequency (ω)
b_33 = inputMatrix(:, 2);      % Second column: b_33 values
b_35 = inputMatrix(:, 3);      % Third column: b_35 values
b_53 = inputMatrix(:, 4);      % Fourth column: b_53 values
b_55 = inputMatrix(:, 5);      % Fifth column: b_55 values

% Initialize τ and preallocate arrays for Ii
tau = 80; 
N = length(frequency); % Number of frequency rows

% Preallocating for integrals
I_33 = zeros(N, 1);
I_35 = zeros(N, 1);
I_53 = zeros(N, 1);
I_55 = zeros(N, 1);

% Calculate integrals (Ii) for each frequency band
for i = 2:N
    % Frequency bounds
    omega_prev = frequency(i-1);
    omega_curr = frequency(i);
    
    % b(ω) values for the current frequency range
    b33_vals = b_33(i-1:i);
    b35_vals = b_35(i-1:i);
    b53_vals = b_53(i-1:i);
    b55_vals = b_55(i-1:i);
    
    % Integrals for each mode (Ii)
    I_33(i) = integral(@(omega) interp1([omega_prev, omega_curr], b33_vals, omega) .* cos(omega * tau), omega_prev, omega_curr);
    I_35(i) = integral(@(omega) interp1([omega_prev, omega_curr], b35_vals, omega) .* cos(omega * tau), omega_prev, omega_curr);
    I_53(i) = integral(@(omega) interp1([omega_prev, omega_curr], b53_vals, omega) .* cos(omega * tau), omega_prev, omega_curr);
    I_55(i) = integral(@(omega) interp1([omega_prev, omega_curr], b55_vals, omega) .* cos(omega * tau), omega_prev, omega_curr);
end

% Calculate damping coefficients B(τ)
B_33 = (2 / pi) * sum(abs(I_33));
B_35 = (2 / pi) * sum(abs(I_35));
B_53 = (2 / pi) * sum(abs(I_53));
B_55 = (2 / pi) * sum(abs(I_55));

% Display the results
fprintf('Damping Coefficients:\n');
fprintf('B_33 = %.4f\n', B_33);
fprintf('B_35 = %.4f\n', B_35);
fprintf('B_53 = %.4f\n', B_53);
fprintf('B_55 = %.4f\n', B_55);
