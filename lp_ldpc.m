function y = lp_ldpc(n,k,d, d_v, d_c)
    %Returns 1 if code is not possible.
    %Valid only for CSS codes where the X and Z stabilizers are identical.
    %Furthermore, the stabilizers have wt d_c
    %qubits connected to d_v stabs of X type and d_v stabilizers of Z type.
    %We assume a dual-containing code of the form Cperp within C.
    %
    %
    %Here V is the weight enumerate of the code Cperp and W is the weight enumerator of the code C
    % V_i is the number of words of weight i in Cperp
    % W_j is the numberof words of weight j in C
    %
    %From the MacWilliams identities, we have
    %       V_i = (1/|C|) \sum_(j=1)^(n) W_j K_i(j).
    %Inequality constraints specified by A.x <= b
    %1.     For all i > d:
    %           sqrt(2^(n+k)) * V_i <= sum(j) V_j K_i(j)
    %           This constraint says that the number of words of weight i in the code C is less than the number of words of weight i in Cperp.
    %
    %2.     Let H parity check matrix be (d_v, d_c) regular.
    %         mu = d_c*(d_v-1) + 1;
    %         L = n*d_v/d_c is number of rows;
    %         For all m <= L/mu,
    %           W_{m * d_c} <= nchoosek(L/mu, m) mu^m
    %         This constraint bounds the number of words in a ball of radius m around a check.
    %
    %         In terms of V, this reads
    %           sum(j) V_j K_{m*d_c}(j) <= sqrt(2^(n+k))*nchoosek(L/mu, m) mu^m
    
    %Equality constraints specified by Aeq.x = beq
    %1. sum(i) V_i = sqrt(2^(n+k))
    %
    %2. V_0 = 1
    % The dual code contains the word 0^n.
    %
    %3. For 1 < i < d: sqrt(2^(n+k)) * V_i = sum(j) V_j K_i(j)
    % This constraint means that for all i less than d,
    % the number of words of weight i in the code Cperp is equal to the number of words of weight i in C.
    
    %Constants:
    L = n*d_v/d_c;
    mu = d_c*(d_v - 1) + 1;
    md = floor(L/mu);
    
    %First n+1 rows of M write W
    M = zeros(n+1+md,n+1);
    for i = 1:n+1
        for j = 1:n+1
            M(i,j) = K(i-1,j-1,n);
        end
    end
    
    %Then subtract W
    for i = 2:n+1
        M(i,i) = M(i,i) - sqrt(2^(n+k));
    end
    
    %Next md rows write LDPC condition
    for m = 1:md
       for j = 1:n+1
            M(n+1+m,j) = -1*K(m*d_c - 1,j-1,n);
        end
    end

    Aeq = zeros(d+1,n+1);
    
    Aeq(1:d,:) = M(1:d,:);
    Aeq(d+1,1) = 1;
        
    beq = zeros(d+1,1);    
    beq(1) = sqrt(2^(n+k));
    beq(d+1) = 1;
    
    A = M(d+1:n+1+md,:);
    b = zeros(n-d+1+md,1);
    for m = 1:md
        b(n-d+1+m) = -1*sqrt(2^(n+k))*bin(L,mu,m);
    end
    
    lb = zeros(n+1,1);
    
    f = zeros(n+1,1);
    
    %options = optimoptions('linprog','Algorithm','dual-simplex','Display','none','OptimalityTolerance',1.0000e-07);
    options = optimoptions('linprog','Algorithm','interior-point','Display','none','OptimalityTolerance',1.0000e-07);
    [s,Z,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,[],options);
    
    y = (-1*exitflag+1)/3;
end