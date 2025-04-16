function y = K(l,i,n)
    %Krawtchuk polynomial over {0,1}^n
    %K(l,i;n) = sum_(a, |a| = l) (-1)^(a.x)
    y = 0;
    for j = 0:l
        y = y + (-1)^(j)*nCk(i,j)*nCk(n-i,l-j);
    end

end

