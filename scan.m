n_max = 25;
n_min = 15;

R = zeros(1,1);
delta = zeros(1,1);
for n = n_min:n_max
    for d = 0:n
        W = lp(n,d);
        Ropt = log2(sum(W))/n;
        delta_opt = d/n;
        R = [R Ropt];
        delta = [delta delta_opt];
    end
end