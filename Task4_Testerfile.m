sizes = [100, 500, 1000];   
k_values = [1, 5, 10];      
num_trials = 3;             

disp("Testing Banded Matrix-Vector Multiplication")

for n = sizes
     disp("Matrix Size: ")
     disp(n)

    for k = k_values
        disp("Band offset: ");
        disp(k);

        A = zeros(n);
        for i = 1:n
            A(i,i) = randn(); 
            if i + k <= n
                A(i, i + k) = randn();  
            end
            if i - k >= 1
                A(i, i - k) = randn();  
            end
        end

        v = randn(n,1);
        v_dense = A * v;
        z = task4(A, v, n);

        % Compare results
        err = norm(v_dense - z(:)) / norm(v_dense);
        disp("Relative Error");
        disp(err);

        % Timing comparison
        t_dense = timeit(@() A * v);
        t_task4 = timeit(@() task4(A, v, n));
        
        disp("Dense Time:")
        disp(t_dense)
        disp("Task 4 Time:")
        disp(t_task4)

    end
end
