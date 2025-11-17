%Task 3 works - Do not mess with it 

function [z] = Task3(A, v, n)

  AA = [];    
  JR = [];    
  JC = [];    
  for i = 1:n
    for j = 1:i    
      if A(i,j) ~= 0
        AA(end+1) = A(i,j);
        JR(end+1) = i;
        JC(end+1) = j;
      end
    end
  end

  %Converts to CSR Form 
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


  z = zeros(n,1);
  for i = 1:n
    k1 = IA(i);
    k2 = IA(i+1) - 1;
   
    for k = k1:k2
      j = JC(k);
      aij = AA(k);
      z(i) = z(i) + aij * v(j);

      if i ~= j
        z(j) = z(j) + aij * v(i);
        
      end
    end
  end

  %disp(IA)
return 

 