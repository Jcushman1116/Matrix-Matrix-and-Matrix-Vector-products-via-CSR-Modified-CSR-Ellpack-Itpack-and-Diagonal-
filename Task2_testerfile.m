sizes = [100, 500,1000];  
density = 0.01;            
elbow_per_row = 5;        

disp("Testing Relaxed CSR Matrix-Vector Multiplication")

for n = sizes
    disp("Matrix Size: ")
    disp(n) 

    
    % Generate random sparse matrix A in COO format
    total_elems = round(n^2 * density);
    JR = randi(n, total_elems, 1);   % Row indices
    JC = randi(n, total_elems, 1);   % Column indices
    AA = randn(total_elems, 1);      % Random nonzero values
    
    % Create full (dense) matrix and vector for validation
    A_dense = full(sparse(JR, JC, AA, n, n));
    v1 = randn(n, 1);
    v2_dense = A_dense * v1;
    v2_relaxed = COO_to_relaxed_CSR_Mv_mult(AA, JR, JC, v1, n, elbow_per_row);
    
    % Compare results
    err = norm(v2_dense - v2_relaxed(:)) / norm(v2_dense);
    disp("Initial product relative error");
    disp(err);

    
    % Add new random elements to A 
    new_elems_per_row = 2;
    added = n * new_elems_per_row;
    JR_add = randi(n, added, 1);
    JC_add = randi(n, added, 1);
    AA_add = randn(added, 1);
    
    JR = [JR; JR_add];
    JC = [JC; JC_add];
    AA = [AA; AA_add];
    
    A_dense = full(sparse(JR, JC, AA, n, n));
    

    % Recompute after insertion
    v2_dense_new = A_dense * v1;
    v2_relaxed_new = COO_to_relaxed_CSR_Mv_mult(AA, JR, JC, v1, n, elbow_per_row);
    
    err_new = norm(v2_dense_new - v2_relaxed_new(:)) / norm(v2_dense_new);
    disp("After insertion Relative Error");
    disp(err_new);

    
    
    % Timing comparison
    t_dense = timeit(@() A_dense * v1);
    t_relaxed = timeit(@() COO_to_relaxed_CSR_Mv_mult(AA, JR, JC, v1, n, elbow_per_row));
    
    disp("Dense Time:")
    disp(t_dense)
    
    disp("Relaxed Time:")
    disp(t_relaxed)
    
end

