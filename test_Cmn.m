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
D = 1; 
L1z = D;
L2z = 0;
L2y = D;
L1 = [0,0,L1z];
L2 = [0,L2y,L2z];
Nl = 10; % Number of unit cells in one direction
N = 2*Nl;
Nn = 2; 

R = 0.1*ones(1,N);


c10 = -1/3*L1;
c20 = 1/3*L1;
c = [];
counter = 1;
for l1 = 1:Nl
    for l2 = 1:1
        c = [c; c10 + l1*L1 + l2*L2; c20 + l1*L1 + l2*L2];
        ind(l1,l2) = counter;
        counter = counter + 1;
    end
end

%%% Plot the geometry
figure, hold on
t = linspace(0,2*pi);
for n = 1:N
    plot(c(n,3)+R(n)*cos(t), c(n,2)+R(n)*sin(t),'k')  
    text(c(n,3),c(n,2),num2str(n))
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
gamma_skin = 0;

matC_static = MakeC_mn(R,c,k0,N_multi);
GCM_static = diag(delta.*v2./vol)*matC_static;

%%% Compute eigenmodes
[evec_static,eval_static] = eig(GCM_static);
[resonances_static,I] = sort(sqrt(diag(eval_static)),'ComparisonMethod','real');
modes_static = evec_static(:,I);

figure
hold on
for j = 1:N
    subplot(5,N/5,j); 
    plot(1:N,real(modes_static(:,j)))
end

%%% Compute the skin capacitance matrix using Multipole

matC_skin = MakeCmn_skin(0,R,c,k0,N_multi);
GCM_skin = diag(delta.*v2./vol)*matC_skin;

%%% Compute eigenmodes
[evec_skin,eval_skin] = eig(GCM_skin);
[resonances_skin,I] = sort(sqrt(diag(eval_skin)),'ComparisonMethod','real');
modes_skin = evec_skin(:,I);

figure
hold on
for j = 1:N
    subplot(5,N/5,j); 
    plot(1:N,real(modes_skin(:,j)))
end

%%% Asymptotic 
matC_asym = makeC_skin_asymp(matC_static,gamma_skin,c,N);
GCM_asym = diag(delta.*v2./vol)*matC_asym;

[evec_asym,eval] = eig(GCM_asym);
modes_asym = evec_asym(:,I);

figure
hold on
for j = 1:N
    subplot(5,N/5,j); 
    plot(1:N,real(modes_asym(:,j)))
end