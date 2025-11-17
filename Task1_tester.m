% % Below are the test problems for Task 1
% % The given inputs are all in-order coordinate form of sparse matrices
% % The vectors to multiply in each case are given in dense non-sparse  form
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Matrix 1
% AA = [10,-1,3,11,5,12,7,13,9,2,3];
% JR = [1,1,2,2,3,3,4,4,5,5,5];
% JC = [1,5,1,2,2,3,3,4,4,1,5]; %Column out of order at the end
% 
% %Vector 1 to Multiply Against Matrix 1
% v = [1; 0; -1; 0; 2];
% n1 = 5; 
% 
% %These should all produce the same results
% 
% disp("Perform the multiplication of Matrix 1 and Vector 1 using the COO to CSR conversion")
% disp(COO_to_CSR_Mv_mult(AA,JR,JC,v,n1))
% 
% disp("Perform the multiplication of Matrix 1 and Vector 1 using the COO to modified CSR conversion")
% disp(COO_to_modified_CSR_Mv_mult(AA,JR,JC,v,n1))
% 
% disp("Perform the multiplication of Matrix 1 and Vector 1 using the COO to modified CSR conversion")
% disp(COO_to_ELL_Mv_mult(AA,JR,JC,v,n1))
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Matrix 2
% AA2 = [1, 5, -2, 4, 1, 2, 1, 2, 3];
% JR2 = [1, 1, 2, 2, 3, 3, 4, 4, 4];
% JC2 = [1, 4, 1, 2, 1, 4, 2, 3, 4];
% 
% Vector 2 to Multiply Against Matrix 2
% v2 = [3; 5; -1; -1];
% n2=4;
% 
% disp("Perform the multiplication of Matrix 2 and Vector 2 using the COO to CSR conversion")
% disp(COO_to_CSR_Mv_mult(AA2,JR2,JC2,v2,n2))
% 
% disp("Perform the multiplication of Matrix 2 and Vector 2 using the COO to modified CSR conversion")
% disp(COO_to_modified_CSR_Mv_mult(AA2,JR2,JC2,v2,n2))
% 
% disp("Perform the multiplication of Matrix 2 and Vector 2 using the COO to modified CSR conversion")
% disp(COO_to_ELL_Mv_mult(AA2,JR2,JC2,v2,n2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Matrix 3
AA3 = [1, 2, 3, 2, -3, 1, 2, 3, -1, 1, 2, 1, -1, 4, 1, 1, -4, 5];
JR3 = [1, 1, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8];
JC3 = [1, 3, 5, 2, 6, 3, 5, 4, 2, 5, 3, 4, 5, 6, 2, 7, 4, 8];

%Vector 3 to Multiply Against Matrix 3
v3 = [-2; 2; -2; 2; -2; 2; -2; 2];
n3=8;

disp("Perform the multiplication of Matrix 3 and Vector 3 using the COO to CSR conversion")
disp(COO_to_CSR_Mv_mult(AA3,JR3,JC3,v3,n3))

disp("Perform the multiplication of Matrix 3 and Vector 3 using the COO to modified CSR conversion")
disp(COO_to_modified_CSR_Mv_mult(AA3,JR3,JC3,v3,n3))

disp("Perform the multiplication of Matrix 3 and Vector 3 using the COO to modified CSR conversion")
disp(COO_to_ELL_Mv_mult(AA3,JR3,JC3,v3,n3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

