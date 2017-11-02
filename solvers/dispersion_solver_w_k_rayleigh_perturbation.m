function [w_d,zeta,PHI_d,t_kloop] = dispersion_solver_w_k_rayleigh_perturbation(kappa,Kf,Cf,Mf,dof_sets,R,n_curves)

% Make Sure Kf and Mf are sparse
density = nnz(Kf)/numel(Kf);
if density<0.2
    Kf = sparse(Kf);
    Mf = sparse(Mf);
end
[~,n_kap] = size(kappa);

%% Form periodicity transformation
% ======================================================================= %
% form Bloch-periodicity matrices using boundary DOF sets
% Output T_per is a structure with fields s0, s1, s2, s3, s12, s13, s23,
% s123. These can be combined to form the Bloch-periodicity transformation
% matrix for any wavevector.
T_per = Periodic_Boundary_Conditions(dof_sets);


fprintf('Commencing Rayleigh Perturbation Solution\n')
fprintf('Model Size: %iDOF (%i after enforcing periodicity)\n',...
    size(T_per.s0,1),size(T_per.s0,2));
%% Dispersion Solution
% ======================================================================= %

n_dof = size(T_per.s0,1);
n_dof_per = size(T_per.s0,2);
if n_curves>n_dof_per
    n_curves = n_dof_per;
end
        
% preallocate Arrays
w_d = zeros(n_curves,n_kap);
zeta = zeros(n_curves,n_kap);
PHI_d = zeros(n_dof,n_curves,n_kap);
t_kloop = zeros(1,n_kap);

% Full FE model Dispersion Solution
tic
for j1 = 1:n_kap

    tstart = toc;

    % wave vector at current k-point
    kvec = kappa(:,j1);
    
    % phase modification at boundaries due to wavevector
    lam  = exp(-1i*R'*kvec);
    lam(end+1:3) = 0;
    T_per_k = T_per.s0 + T_per.s1*lam(1) + T_per.s2*lam(2) + T_per.s3*lam(3)...
            + T_per.s12*lam(1)*lam(2) + T_per.s23*lam(2)*lam(3) + T_per.s13*lam(1)*lam(3)...
            + T_per.s123*lam(1)*lam(2)*lam(3);
        
    % apply periodicity transformation to mass and stiffness matrices
    K_hat = (T_per_k'*Kf*T_per_k);
    M_hat = (T_per_k'*Mf*T_per_k);
    C_hat = (T_per_k'*Cf*T_per_k);

    %symmetrize matrices
    K_hat = (1/2)*(K_hat+K_hat');
    M_hat = (1/2)*(M_hat+M_hat');
    C_hat = (1/2)*(C_hat+C_hat');
    
    % Compute Undamped eigensolution 
    if n_dof_per/n_curves<=5
%     if n_dof_per<=800
        [PHI_und,L_und] = eig(full(K_hat),full(M_hat)); 
    else

        % eigensolution (if zero frequency modes exist, provide a small
        % perturbation shift and try again)
        try
            [PHI_und,L_und] = eigs(K_hat,M_hat,n_curves,'sm');
        catch
            [PHI_und,L_und] = eigs(K_hat,M_hat,n_curves,0.01);
        end
    end
    w_n = sqrt(sort(abs(diag(L_und))));
    
    % Compute damped eigensolulution using rayleigh perturbation
    [PHI_temp,D_damped]=RPM(M_hat,C_hat,K_hat,PHI_und,L_und,[],[]);
    PHI_d(:,:,j1) = T_per_k*PHI_temp;
    
    zeta(:,j1) = imag(D_damped)./w_n;
    w_d(:,j1)=w_n.*sqrt(1-(zeta(:,j1).^2));

    t_kloop(j1) = toc-tstart;
    
    % display loop timing info
    fprintf('k-point %i of %i, solution time: %4.2f\n',j1,n_kap,t_kloop(j1))
end

w_d = real(w_d);
zeta = real(zeta);