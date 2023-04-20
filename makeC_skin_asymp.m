function C = makeC_skin_asymp(C1,gamma_skin,c)
C = C1;
cz = c(:,3);
N = length(cz);
%     for l = 1:N_skin
%         for j = 1:N_skin
%             C(l,j)=C1(l,j)*(1-eps_ski*(cz(l)-cz(j)));
%         end
%     end
    for l = 1:N
        for j =1:N
             C(l,j)=C1(l,j)*(1+gamma_skin*(cz(l)-cz(j)));
        end
    end
end