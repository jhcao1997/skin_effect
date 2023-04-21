%
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% Define parameters
%%% Geometry
% Define the length of the unit cell
N = 25;
cx = [0:1:N-1]';
cy = zeros(N,1);
cz = zeros(N,1);
c = [cx cy cz];

R = 0.1*ones(N,1)';
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
delta = 1./9000;
v2 = ones(1,N);
N_multi = 2;

%%% Compute the static capacitance matrix
% Maximum order for multipole expansion (n = 0, 1, ..., N_multi)
% If we use higher order, then accuracy improves. Usually 0 is sufficiently large.
gamma_skin = 0.1;

matC_static = MakeC_mn(R,c,k0,N_multi);
GCM_static = diag(delta.*v2./vol)*matC_static;

%%% Compute eigenmodes
[evec_static,eval_static] = eig(GCM_static);
[resonances_static,I] = sort(sqrt(diag(eval_static)),'ComparisonMethod','real');
modes_static = evec_static(:,I);

figure
title("modes_static")
hold on
for j = 1:N
    subplot(5,N/5,j); 
    plot(1:N,real(modes_static(:,j)))
end

%%% Compute the skin capacitance matrix using Multipole

matC_skin = MakeCmn_skin(gamma_skin,R,c,k0,N_multi);
GCM_skin = diag(delta.*v2./vol)*matC_skin;

%%% Compute eigenmodes
[evec_skin,eval_skin,eigen_left] = eig(GCM_skin);
[resonances_skin,I] = sort(sqrt(diag(eval_skin)),'ComparisonMethod','real');
modes_skin = evec_skin(:,I);

figure
title("modes_skin")
hold on
for j = 1:N
    subplot(5,N/5,j); 
    hold on
    plot(1:N,real(exp(-gamma_skin*cx).*modes_skin(:,j)))
%     plot(1:N,real(eigen_left(:,j)))

%     plot(1:N,imag(modes_skin(:,j)))
end

%%% Asymptotic 
matC_asym = makeC_skin_asymp(matC_static,gamma_skin,c);
GCM_asym = diag(delta.*v2./vol)*matC_asym;

[evec_asym,eval] = eig(GCM_asym);
modes_asym = evec_asym(:,I);
% 
% figure
% hold on
% for j = 1:N
%     subplot(5,N/5,j); 
%     plot(1:N,real(modes_asym(:,j)))
% end


% compute the degree of localization
N_deg = 15;
gammas = linspace(0,0.8,N_deg);
average_deg = zeros(N_deg,1);
max_deg = zeros(N_deg,1);
min_deg= zeros(N_deg,1);
for i = 1:N_deg
    degs = zeros(N,1);
    matC_skin = MakeCmn_skin(gammas(i),R,c,k0,N_multi);
    GCM_skin = diag(delta.*v2./vol)*matC_skin;
    [evec_deg,~] = eig(GCM_skin);
    for j = 1:N
        degs(j) = norm(evec_deg(:,j),Inf)/norm(evec_deg(:,j),2);
    end
    average_deg(i) = mean(degs);
    max_deg(i) = max(degs);
    min_deg(i) = min(degs);
end
figure
hold on
plot(gammas,average_deg);
plot(gammas,max_deg);
plot(gammas,min_deg);
legend('average','max','min')
xlabel('\gamma')
ylabel('Degrees of localization')
%%% Compute eigenmodes
[evec_skin,eval_skin] = eig(GCM_skin);
[resonances_skin,I] = sort(sqrt(diag(eval_skin)),'ComparisonMethod','real');
modes_skin = evec_skin(:,I);

