function y = lp(n,k,d)
    %Valid only for CSS codes where the X and Z stabilizers are identical.
    %Furthermore, the stabilizers have wt d_c
    %qubits connected to d_v stabs of X type and d_v stabilizers of Z type.
    %
    %
    %Inequality constraints specified by A.x ? b
    %1.     For all i ? d:
    %           sqrt(2^(n+k)) * W'_i ? sum(j) W'_j K_i(j)
    %
    %2.     Let H parity check matrix be (d_v, d_c) regular.
    %         µ = d_c(d_v-1) + 1;
    %         L = n*d_v/d_c is number of rows;
    %         For all m ? L/µ,
    %           W_{md_c} ? nchoosek(L/µ, m) µ^m
    %
    %         In terms of the W', this reads
    %           sum(j) W'_j K_{md_c}(j) <= sqrt(2^(n+k))*nchoosek(L/µ, m) µ^m
    
    %Equality constraints specified by Aeq.x = beq
    %1. sum(i) W'_i = sqrt(2^(n+k))
    %
    %2. W'_0 = 1
    %
    %3. For 1 ? i ? d: sqrt(2^(n+k)) * W'_i = sum(j) W'_j K_i(j) 
    
    %First n+1 rows of M write W
    M = zeros(n+1,n+1);
    for i = 1:n+1
        for j = 1:n+1
            M(i,j) = K(i-1,j-1,n);
        end
    end
    
    %Then subtract W
    for i = 2:n+1
        M(i,i) = M(i,i) - sqrt(2^(n+k));
    end

    Aeq = zeros(d+1,n+1);
    
    Aeq(1:d,:) = M(1:d,:);
    Aeq(d+1,1) = 1;
        
    beq = zeros(d+1,1);    
    beq(1) = sqrt(2^(n+k));
    beq(d+1) = 1;
    
    A = M(d+1:n+1,:);
    b = zeros(n-d+1,1);
    
    lb = zeros(n+1,1);
    
    f = zeros(n+1,1);
    
    %options = optimoptions('linprog','Algorithm','dual-simplex','Display','none','OptimalityTolerance',1.0000e-07);
    options = optimoptions('linprog','Algorithm','interior-point','Display','none','OptimalityTolerance',1.0000e-07);
    [s,Z,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,[],options);
    
    y = (-1*exitflag+1)/3;
end