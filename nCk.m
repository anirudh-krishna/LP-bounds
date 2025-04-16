function y = nCk(n,k)
    if k > n
        y = 0;
    else
        y = exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1));
    end
end
