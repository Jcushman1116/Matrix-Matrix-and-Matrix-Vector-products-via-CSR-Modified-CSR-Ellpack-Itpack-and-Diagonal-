function [z] = task4(A,v,n)
    
    Diag = zeros(n, 3);    
   
    k = 0;  
    for i = 1:n-1
        if A(i, i+1) ~= 0   
            k = i;            
            break;
        elseif A(i+1, i) ~= 0  
            k = i;            
            break;
        end
    end

    IOFF = [-k, 0, k];

    for i = 1:n
        Diag(i, 2) = A(i, i);
    end

    for i = (k+1):n
        Diag(i, 1) = A(i, i-k);
    end

    for i = 1:(n-k)
        Diag(i, 3) = A(i, i+k);
    end

    z = zeros(n,1);
    for i =1:n 
        for j = 1:length(IOFF)
            if i + IOFF(j) >= 1 && i + IOFF(j) <= n
                z(i) = z(i) + Diag(i,j) * v(i + IOFF(j));
            end
        end
    end 

return

   