function C = makeC_skin(C1,eps_ski,c,N_skin)
C = C1;
N = length(c);
    for l = 1:N_skin
        for j = 1:N_skin
            C(l,j)=C1(l,j)*(1+0.5i*eps_ski*(c(l)-c(j)));
        end
    end
%     for l = N_skin+1:N
%         for j =1:N_skin
%              C(l,j)=C1(l,j)*(1+0.5i*eps_ski*(c(l)-c(j)));
%         end
%     end
%     for l = 1:N_skin
%         for j =N_skin+1:N
%              C(l,j)=C1(l,j)*(1+0.5i*eps_ski*(c(l)-c(j)));
%         end
%     end
    C_test = C-C1;
end