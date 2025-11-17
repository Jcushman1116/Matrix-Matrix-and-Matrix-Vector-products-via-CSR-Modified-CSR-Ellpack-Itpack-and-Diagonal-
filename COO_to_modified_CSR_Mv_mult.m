%Task 1 part 2

function [z] = COO_to_modified_CSR_Mv_mult(AA, JR, JC, v, n)

diagonal_elements = zeros(1, n);
nondiagonal_elements = [];
nondiagonal_elements_columns = [];
rowpointers = zeros(1, n + 1);
z = zeros(n, 1);

% input diagonal and non-diagonal elements (and associated columns)
for i = 1:length(AA)
    if JR(i) == JC(i)
        diagonal_elements(JR(i)) = AA(i);
    else 
        nondiagonal_elements = [nondiagonal_elements, AA(i)];
        nondiagonal_elements_columns = [nondiagonal_elements_columns, JC(i)];
    end 
end 

%find row pointers and fill array 
p = 1; 
for i = 1:n
    rowpointers(i) = p; 
    p = p + sum(JR ~= JC & JR == i); % Count non-diagonal elements for row i
end 
rowpointers(n + 1) = p; 

% Append arrays to form modified CSR storage
modified_AA = [diagonal_elements, 0, nondiagonal_elements];
JA = [rowpointers, nondiagonal_elements_columns]; 

disp(modified_AA)
disp(JA)

% Matrix Vector Mult
for i = 1:n 
    z(i) = modified_AA(i)*v(i);
    for j = JA(i):JA(i+1)-1
        col_index = JA(n+1+j);
        z(i) = z(i) + modified_AA(j+n+1)*v(col_index);
    end
end

end 
