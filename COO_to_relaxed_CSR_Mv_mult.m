

function [z] = COO_to_relaxed_CSR_Mv_mult(AA,JR,JC, v, n, elbow)

    sorted_JR = sort(JR);
    IA = zeros(n+1,1);
    IA(1) = 1;

    for i = 1:n
        rowsum = sum(sorted_JR == i);
        IA(i+1) = IA(i) + rowsum;
    end

    AA_CSR = zeros(length(AA),1);
    JC_CSR = zeros(length(JC),1);
    p = 1;
    for i = 1:n
        for j = 1:length(JR)
            if JR(j) == i
                AA_CSR(p) = AA(j);
                JC_CSR(p) = JC(j);
                p = p + 1;
            end
        end
    end

    elbowroom = n * elbow;
    AA_relaxed = zeros(length(AA_CSR) + elbowroom, 1);
    JC_relaxed = zeros(length(JC_CSR) + elbowroom, 1);
    IA_relaxed = zeros(n+1,1);
    IA_relaxed(1) = 1;

    h = 1;
    for i = 1:n
        k1 = IA(i);
        k2 = IA(i+1)-1;
        rowlength = k2 - k1 + 1;

        AA_relaxed(h:h+rowlength-1) = AA_CSR(k1:k2);
        JC_relaxed(h:h+rowlength-1) = JC_CSR(k1:k2);

        h = h + rowlength + elbow;
        IA_relaxed(i+1) = h;
    end
    IA = IA_relaxed; 

    for i = 1:n
        k1 = IA(i);
        k2 = IA(i+1)-elbow-1;
        if k1 <= k2
            z(i) = dot(AA_relaxed(k1:k2), v(JC_relaxed(k1:k2)));
        end
    end
end
