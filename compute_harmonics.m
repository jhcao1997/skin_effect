function f_out = compute_harmonics(fun,N_multi)
    f_out = zeros(lin_ind(N_multi,N_multi),1);
    ind = 1;
    for n = 0:N_multi
        for m = -n:n
            integrand = @(theta, phi) fun(theta,phi)*sin(theta)*conj(harmonicY(n,m,theta,phi));
            f_out(ind) = int_trapez(integrand,200,0,pi,0,2*pi);
            ind = ind+1;
        end
    end
end

function res = int_trapez(integrand,NN,xmin,xmax,ymin,ymax)
    xs = linspace(xmin,xmax,NN);
    ys = linspace(ymin,ymax,NN);
    mesh = zeros(NN,NN);
    for i = 1:NN
        for j =1: NN
            mesh(i,j) = integrand(xs(i),ys(j));
        end
    end
    res = trapz(ys,trapz(xs,mesh,2));
end