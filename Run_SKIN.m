
% Number of Resonators
k = 20;
ll = 1;
N_skin = k*ll;
N = N_skin+3;

% Perturbed case, change eps = 0 to obatian the unperturbed case
eps = 0;
Omega = 0.2;
eps_skin = 1i*0.2;
[V,f_exp,f_exp_im,ind,c,M_0] = produce_edge_supercell(k,ll,eps,Omega,eps_skin,N);

% Plot all the bloch modes
figure 
hold on 
for j = 1:N
    W = exp(-eps_skin*1i*c(:,1)).*V(1:N,j);
%    W = exp(0*1i*c(:,1)).*V(1:N,j);
     subplot(5,ceil(N/5),j)
     hold on
    plot(1:N,real(W),'b');
%     plot(1:N,imag(W),'r');
    xline(N_skin)
%     legend('Re','Im')
%     xlabel(strcat('The serial number of the resonators',' ,eps=',num2str(eps)))
end

function [V,f_exp,f_exp_im,ind,c,M_0] = produce_edge_supercell(k,ll,eps,Omega,eps_skin,N)
    N_skin = k*ll;
    % radii
    R = 0.1*ones(1,N);

    %%% Material parameters
    rho0 = 1e3;             % background material
    kappa0 = 1e3;           % background material
    v = sqrt(kappa0/rho0);  % speed of sound in water

    rho_b = 1;            % density of resonators  
    kappa_b = 1;          % bulk modulus of resonators
    v_b = sqrt(kappa_b/rho_b);  % speed of sound resonators

    
    % supercell
    D = 1; 
    L1x = D*sqrt(3);
    L2x = D*sqrt(3)/2;
    L2y = D*3/2;
    L1 = [L1x,0];
    L2 = [L2x,L2y];

    eps = 3*R(1); %eps = 3*R;

    L_unit = sqrt(3);
    cx = [0];


    for ii = 1: N-1
        cx = [cx ; ii*L_unit];
    end

    phase_shift= zeros(N,1);   

    cy = zeros(N,1);
    c = [cx cy];
    % Plot the geometry
    figure, hold on
    t = linspace(0,2*pi);
    for n = 1:N
        plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
        text (cx(n), cy(n), num2str(n))
    end
    
    daspect([1 1 1])
    hold off
    close

    rhot = @(t) rho_b./(1 + eps*cos(Omega*t + phase_shift)); 

    sqrtkappat = @(t) sqrt(kappa_b)*ones(N,1); % ./[1 + eps*sin(t + pi*(1:N))];
    dkappa = @(t) 0*R;
    d2kappa = @(t) 0*R;
    w3 = @(t) 0*(1:N);

    T = 2*pi/Omega;
    steps = 1000;

    % High contrast parameters \delta
    delta=rho_b/rho0;

    %% Compute the resonances using the capacitance matrix, which is approximated using the multipole expansion method
    C = capacitance_2D(c,R,rho0,rho_b,kappa0,kappa_b,delta);
    C = makeC_skin(C,eps_skin,cx,N_skin);
    vol = 4/3*pi*R.^3;
    vol = vol';
    CoeffMat = @(t) makeM(t,delta,kappa0,rho0,vol,C,rhot, sqrtkappat, w3);
    M_0 = delta*kappa0/rho0*inv(diag(vol))*C;
    %% Solve for Psi
    [Tspan, X_fundamental] = HillSolver(CoeffMat,T,steps); % X_fundamental = [Psi; dPsi/dt]
    [TOUT, X_bloch,V_mode,V,D] = Hill2BlochModes(CoeffMat,T,steps);
    
    
    %% Solve Floquet exp
    X_T = squeeze(X_fundamental(end,:,:));
    [V,d] = eig(X_T,'vector');
    [f_exp,ind] = sort(imag(log(d)/T),'descend');
    f_exp_im = real(log(d(ind))/T);
    V = V(:,ind);
end