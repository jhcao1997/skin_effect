%
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
%% Define parameters
%%% Geometry
% Define the length of the unit cell
N = 10;
cx = [0:1:N-1]';
cy = zeros(N,1);
cz = zeros(N,1);
c = [cx cy cz];

R = 0.3*ones(N,1)';
%%% Plot the geometry
figure, hold on
t = linspace(0,2*pi);
for n = 1:N
    plot(c(n,1)+R(n)*cos(t), c(n,2)+R(n)*sin(t),'k')  
    text(c(n,1),c(n,2),num2str(n))
end
daspect([1 1 1])
hold off
close 

vol = 4*pi*R.^3/3;
k0 = 0.0000001;
delta = 10^(-5);
v2 = ones(1,N);
N_multi = 2;
%%% Compute the static capacitance matrix
% Maximum order for multipole expansion (n = 0, 1, ..., N_multi)
% If we use higher order, then accuracy improves. Usually 0 is sufficiently large.
gamma_skin = 1;

matC_static = MakeC_mn(R,c,k0,N_multi);
GCM_static = diag(delta.*v2./vol)*matC_static;

%%% Compute eigenmodes
[evec_static,eval_static] = eig(GCM_static);
[resonances_static,I] = sort(sqrt(diag(eval_static)),'ComparisonMethod','real');
modes_static = evec_static(:,I);


% compute the normalization factor
fun = @(theta,phi,r) exp(gamma_skin*r*sin(theta)*cos(phi))*r^2*sin(theta);
int_A = int_trapez_3(fun,200,0,pi,0,2*pi,0,R(1));
A_norm = exp(gamma_skin*cx)*int_A;

%%% Compute the skin capacitance matrix using Multipole

matC_skin = MakeCmn_skin(gamma_skin,R,c,k0,N_multi);

GCM_skin = diag(delta.*v2./A_norm')*matC_skin;


%% compute winding number
f = symbol(GCM_skin,14);
thetas = linspace(0,2*pi,100);
fs = [];
for theta = thetas
    fs = [fs f(exp(1i*theta))];
%     scatter(real(fs),imag(fs),[],linspace(0,2*pi,length(fs)))
end
figure
H = arrowPlot(real(fs), imag(fs), 'number', 10,'color','k');
hold on
plot(real(eig(GCM_skin)),imag(eig(GCM_skin)),'*','Color','k')

k = 10;
GCM_k = zeros(size(GCM_skin));
for i = -k:k
    GCM_k = GCM_k + diag(diag(GCM_skin,i),i);
end

k = floor((N-1)/2);
means = zeros(k,N);
GCM_k = diag(diag(GCM_skin));
for i = 1:10
    GCM_k = GCM_k + diag(diag(GCM_skin,i),i);
    GCM_k = GCM_k + diag(diag(GCM_skin,-i),-i);
    [evec_skin,eval_skin,eigen_left] = eig(GCM_k);
    means(i,:) = mean(abs(evec_skin),2);
end

errors= zeros(1,k);
for i= 1:k
    errors(i) = norm(means(i,:) -means(k,:) )/norm(means(k,:));
end

figure
plot([1:k],errors,'k')

figure 
hold on
for i = 10:k
    plot(means(i,:))
end
legendStrings = "k = " + string([10:k]);
legend(legendStrings)

%%% Compute eigenmodes
[evec_skin,eval_skin,eigen_left] = eig(GCM_k);
[evec_skin2,eval_skin2,eigen_left2] = eig(GCM_skin);

[resonances_skin,I] = sort(sqrt(diag(eval_skin)),'ComparisonMethod','real');
[resonances_skin,I2] = sort(sqrt(diag(eval_skin2)),'ComparisonMethod','real');

modes_skin = evec_skin(:,I);
modes_skin2 = evec_skin2(:,I);

figure
for j = 1:20
    subplot(ceil(20/5),5,j)
    hold on
    plot(1:N,real(modes_skin(:,j)),'k');
    plot(1:N,real(modes_skin2(:,j)),'.-k')
end


% prediction of the mathematical model 
% modes_math = zeros(N,N);
% for k = 1:N
%     for i = 1:N
%     modes_math(i,k) = (GCM_skin(2,1)/GCM_skin(1,2))^(i/2)*sin(i*k*pi/(N+1));
%     end
%     modes_math(:,k) = modes_math(:,k)/norm(modes_math(:,k));
% end



figure
title("Eigenmodes for \gamma =" + num2str(gamma_skin))
hold on
%         plot(cx,mean(abs(modes_math),2),'k','linewidth',3)
        plot(cx,mean(abs(modes_skin2),2),'r','linewidth',3)
%         plot(cx,mean(abs(modes_skin),2),'g','linewidth',3)
for j = 1:N
%     plot(1:N,real(modes_skin(:,j)),'b')
plot(cx,real(modes_skin2(:,N-j+1)),'color', [.5 .5 .5])

% %     plot(1:N,real(eigen_left(:,j)))
% %     plot(1:N,imag(modes_skin(:,j)))
end
        plot(cx,mean(abs(modes_skin2),2),'r','linewidth',3)

legend('Average of the absolute amplitudes','Eigenmodes')


