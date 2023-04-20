f = @(x,y) fun_g(x,y);

integral2(f,0,1,0,1);
function h = fun_g(x,y)
    h = x.^2+y/3;
end