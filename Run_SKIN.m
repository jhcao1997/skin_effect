% close all
% v = VideoWriter('skinvideo.mp4');
% open(v);
% for i = 0:10
%     gamma = -0.5*i;

N = 50;
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
N_multi = 3;
%%% Compute the static capacitance matrix
% Maximum order for multipole expansion (n = 0, 1, ..., N_multi)
% If we use higher order, then accuracy improves. Usually 0 is sufficiently large.
gamma_skin = -1;

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

%%% Plot the matrix entry to verify Lemma 5.5
GCM_skin = diag(delta.*v2./A_norm')*matC_skin;

figure
hold on
for i = 1:N
    for j =1:N
        scatter((abs(i-j)),abs(GCM_skin(i,j)),'*','black')
    end
end

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Distance between the indices $|i-j|$','interpreter','latex','FontSize',17)
ylabel('Value of the coefficients $|C_{ij}|$','interpreter','latex','FontSize',17)
set(gca,'ticklabelinterpreter','latex')
set(gca, 'FontSize',17)

%% compute winding number

% 
f = symbol(GCM_skin,10);
thetas = linspace(0,2*pi,100);
fs = [];
for theta = thetas
    fs = [fs f(exp(1i*theta))];
%     scatter(real(fs),imag(fs),[],linspace(0,2*pi,length(fs)))
end
ev = eig(GCM_skin);

figure
hold on
H = arrowPlot(real(fs), imag(fs), 'number', 10,'color','k');
hold on
plot(real(eig(GCM_skin)),imag(eig(GCM_skin)),'*','Color','k')
plot(real(ev(1)),imag(ev(1)),'*','Color','r')
xlabel('Real part','FontSize',18,'Interpreter','latex')
ylabel('Imaginary part','FontSize',18,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',23)

GCM_k = zeros(N);
k = 10;
for i = -k:k
    GCM_k = GCM_k + diag(diag(GCM_skin,i),i);
end


%%% Compute eigenmodes
[evec_skin,eval_skin,eigen_left] = eig(GCM_k);
[evec_skin2,eval_skin2,eigen_left2] = eig(GCM_skin);

[resonances_skin,I] = sort(sqrt(diag(eval_skin)),'ComparisonMethod','real');
[resonances_skin,I2] = sort(sqrt(diag(eval_skin2)),'ComparisonMethod','real');

modes_skin = evec_skin(:,I);
modes_skin2 = evec_skin2(:,I2);

% eval_skin = eval_skin(:,I);
% eval_skin2 = eval_skin(:,I2);

figure
for j = 1:20
    subplot(ceil(20/5),5,j)
    hold on
%     plot(1:N,real(modes_skin(:,j)),'b','LineWidth',0.25)
    plot(1:N,real(modes_skin2(:,j)),'k','LineWidth',0.25)
    set(gca,'TickLabelInterpreter','latex','FontSize',10)
    if j == 18
        xlabel('Position of the resonators','Interpreter','Latex')
    end
end
% 

% prediction of the mathematical model 
% modes_math = zeros(N,N);
% for k = 1:N
%     for i = 1:N
%     modes_math(i,k) = (GCM_skin(2,1)/GCM_skin(1,2))^(i/2)*sin(i*k*pi/(N+1));
%     end
%     modes_math(:,k) = modes_math(:,k)/norm(modes_math(:,k));
% end



fig=figure;
hold on
%         plot(cx,mean(abs(modes_math),2),'k','linewidth',3)
%         plot(cx,mean(abs(modes_skin),2),'g','linewidth',3)
for j = 1:N
    plot(1:N,real(modes_skin2(:,j)),'k')
end
% plot(cx,mean(abs(modes_skin2),2),'r','linewidth',3)
% plot(cx,real(modes_skin2(:,1)),'color', 'r')

xlabel('Position of the resonators','FontSize',25)
ylabel('')
set(gca,'TickLabelInterpreter','latex','FontSize',25)
frame = getframe(fig);
writeVideo(v,frame);
% end




