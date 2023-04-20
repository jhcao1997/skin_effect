function f_out = compute_harmonics(fun,Ri,N_multi)
    f_out = zeros(lin_ind(N_multi,N_multi),1);
    ind = 1;
    for n = 0:N_multi
        for m = -n:n
            integrand = @(theta, phi) Ri^2*fun(theta,phi).*sin(theta).*conj(harmonicY(n,m,theta,phi));
            f_out(ind) = trapez(integrand,0,pi,0,2*pi);
            ind = ind+1;
        end
    end
end

function res = trapez(fun,xmin,xmax,ymin,ymax)
    N = 200;
    x= linspace(xmin,xmax,N);
    y = linspace(ymin,ymax,N);
    xmid = (x(1:end-1)+x(2:end))/2;
    ymid = (y(1:end-1)+y(2:end))/2;
    res = 0;
    for i = 1:length(xmid)
        for j =1:length(ymid)
            res = res + fun(xmid(i),ymid(j));
        end
    end
    res = res * (xmax-xmin)*(ymax-ymin)/ (N^2);
end