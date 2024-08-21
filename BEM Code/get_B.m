function B_j = get_B(r, n, area, r_mod)
    % This function calculates the value of B for each point j in the domain
    % Inputs:
    %   r:      n x n x 3 matrix (jth point, ith reference, [x, y, z] components)
    %   n:      n x n x 3 matrix (normal vector components [x, y, z])
    %   area:   n-element vector containing areas of n panels
    %   r_mod:  n x n matrix containing the mod values of r vectors
    % Outputs:
    %   B_j:    n-element vector containing B values for each point j
    
    % Get the number of points/panels
    [n_points, n_refs, ~] = size(r);
    
    % Initialize the output vector B_j
    B_j = zeros(n_points, 1);
    
    % Loop over each reference point and each domain point
    for j = 1:n_points
        b_i_sum = 0; % Initialize sum for each j
        for i = 1:n_refs
            % Skip the case where i == j
            if i == j
                b_i = 0;
            else
                % Extract r_j_i vector (r(j,i,:)) and n_i vector (n(j,i,:))
                r_j_i = squeeze(r(j,i,:));  % [r_x, r_y, r_z]
                n_i = squeeze(n(j,i,:));    % [n_x, n_y, n_z]
                
                % Calculate the dot product of r_j_i and n_i
                dotProduct = dot(r_j_i, n_i);
                
                % Calculate b_i
                b_i = (dotProduct / r_mod(j,i)^2) * area(i);
            end
            
            % Accumulate the sum for B_j
            b_i_sum = b_i_sum + b_i;
        end
        
        % Assign the computed sum to B_j
        B_j(j) = b_i_sum;
    end
end

% Example dimensions
%n = 3;

% Example r matrix (n x n x 3)
r = rand(20, 20, 3);

% Example n matrix (n x n x 3)
n = rand(20, 20, 3);

% Example area vector
area = rand(1, 20);

% Example r_mod matrix (n x n)
r_mod = sqrt(sum(r.^2, 3));  % This is one way to calculate the modulus

% Calculate B_j
B_j = get_B(r, n, area, r_mod);

% Display the result
disp(B_j);
