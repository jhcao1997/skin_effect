function C = makeC_skin_asymp(C1,gamma_skin,c)
C = C1;
cx = c(:,1);
N = length(cx);
    for l = 1:N
        for j =1:N
             C(l,j)=C1(l,j)*(exp(gamma_skin*(cx(l)-cx(j))));
        end
    end
end