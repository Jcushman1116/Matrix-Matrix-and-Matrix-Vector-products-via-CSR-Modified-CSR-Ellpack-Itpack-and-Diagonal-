
sizes = [100, 500, 1000];    
density = 0.01;             
num_trials = 3;             

disp("Testing Symmetric CSR Matrix-Vector Multiplication")

for n = sizes
    disp("Matrix Size: ")
    disp(n)

    % Generate random sparse symmetric matrix (lower + diag)
    total_elems = round((n^2 * density) / 2);   % half since symmetric
    JR = randi(n, total_elems, 1);   % row indices
    JC = randi(n, total_elems, 1);   % col indices
    AA = randn(total_elems, 1);      % random nonzeros

    mask = JR < JC;
    JR(mask) = JC(mask);

    A_lower = sparse(JR, JC, AA, n, n);
    A_sym = A_lower + triu(A_lower.',1); % reflect upper
    A_dense = full(A_sym);

    v1 = randn(n,1);
    v2_dense = A_dense * v1;
    v2_sym = Task3(A_dense, v1, n);

    % Compare results
    err = norm(v2_dense - v2_sym(:)) / norm(v2_dense);
    disp("Relative Error");
    disp(err);

    t_dense = timeit(@() A_dense * v1);
    t_sym   = timeit(@() Task3(A_dense, v1, n));

    disp("Dense Time:")
    disp(t_dense)
    
    disp("Symmetric Time:")
    disp(t_sym)

end

