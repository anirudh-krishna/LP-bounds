function y = bin(L,mu,m)
    prod = 1;
    for i = 1:m
        prod = prod * ((L/mu) - (i-1))* mu/(m-(i-1));
    end
    y = prod;
end