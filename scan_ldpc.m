n_max = 500;
n_min = 100;

k = 1;
interval = 5

N = zeros(1,1);
K = zeros(1,1);
dist = zeros(1,1);

%Search for existence of code between n_min and n_max in intervals of interval
d_v = 2;
d_c = 5;
for n = n_min:interval:n_max
    if mod(n,5) == 0
        disp(n)
    end
    for k = 1:n
        d = floor(0.01*n);
        feasible = lp_ldpc(n,k,d,d_v,d_c);
        if feasible == 1
            N = [N n];
            K = [K k];
            break
        end
    end
end
