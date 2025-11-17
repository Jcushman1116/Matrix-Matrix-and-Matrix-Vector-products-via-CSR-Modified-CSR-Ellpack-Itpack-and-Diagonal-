%Task 1 part 1 

function[z] = COO_to_CSR_Mv_mult(AA,JR,JC,v,n)

  z=zeros(n,1);  
 
  %Gather elements of IA
  IA(1)=1; 
  for i=1:n
       rowsum = 0; 
       for j=1:length(JR)
           if JR(j) == i
               rowsum = rowsum+1;
           end 
       end
       IA(i+1) = IA(i)+rowsum;
  end 
 
  %preforms multiplication of Mv using dot product method from notes
  for i = 1:n 
      k_1 = IA(i);
      k_2 = IA(i+1)-1;
      z(i) = dot(AA(k_1:k_2),v(JC(k_1:k_2)));
  end 

  disp(AA)
  disp(JC)
  disp(IA)

return 

%Works during testing