A = [1 1; -1 -1];
b = [1 -5];

f = [-1 -1];

options = optimoptions('linprog','Algorithm','dual-simplex','Display','none','OptimalityTolerance',1.0000e-07);
[s,Z,exitflag,output,lambda] = linprog(f,A,b,[],[],[],[],options);