%Task 1 part 3 

function [z] = COO_to_ELL_Mv_mult(AA,JR,JC,v,n)

z=zeros(n,1);
rowcounter = zeros(n,1); 
rowpointer = zeros(n,1);

% Determine # nonzero elements in A
for i =1:length(JR)
    rowcounter(JR(i)) = rowcounter(JR(i)) + 1; 
end 

num_nonzeros = max(rowcounter); 
COEF = zeros(n,num_nonzeros); 
JCOEF = zeros(n, num_nonzeros);

%Determine elements of COEF & JCOEF
for k=1:length(AA)
    i = JR(k);
    j = rowpointer(i) +1; 
    COEF(i,j) = AA(k);
    JCOEF(i,j) = JC(k);
    rowpointer(i) = j; 
end

%Matrix Vector mult
for i = 1:n 
    for j = 1:num_nonzeros
        if JCOEF(i,j) > 0
            z(i) = z(i) + COEF(i,j)*v(JCOEF(i,j));
        end
    end
end

disp(COEF)
disp(JCOEF)
end
