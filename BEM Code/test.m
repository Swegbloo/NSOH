% Example dimensions
n = 3;

% Example r matrix (n x n x 3)
r = rand(10, 10, 3);

% Example n matrix (n x n x 3)
n = rand(10, 10, 3);

% Example area vector
area = rand(1, 10);

% Example r_mod matrix (n x n)
r_mod = sqrt(r(:,:,1).^2+r(:,:,2)^2+r(:,:,3)^2);  % This is one way to calculate the modulus

% Calculate A_j_i
A_j_i = get_inf_coef(r, n, area, r_mod);

% Display the result
disp(A_j_i);