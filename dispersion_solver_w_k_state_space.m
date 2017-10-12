function [lambda,PHI,t_kloop] = dispersion_solver_w_k_state_space(kappa,Kf,Cf,Mf,dof_sets,R,n_curves,shift)

% Make Sure Kf and Mf are sparse
density = nnz(Kf)/numel(Kf);
if density<0.2
    Kf = sparse(Kf);
    Mf = sparse(Mf);
    Cf = sparse(Cf);
end
[~,n_kap] = size(kappa);

%% Form periodicity transformation
% ======================================================================= %
% form Bloch-periodicity matrices using boundary DOF sets
% Output T_per is a structure with fields s0, s1, s2, s3, s12, s13, s23,
% s123. These can be combined to form the Bloch-periodicity transformation
% matrix for any wavevector.
T_per = Periodic_Boundary_Conditions(dof_sets);

%% Dispersion Solution
% ======================================================================= %

n_dof = size(T_per.s0,1);
n_dof_per = size(T_per.s0,2);
        
% preallocate Arrays
% lambda = zeros(n_curves,n_kap);
% PHI = zeros(n_dof_per,n_curves,n_kap);
t_kloop = zeros(1,n_kap);
% lambda_full = nan(2*n_dof_per,n_kap);

% Full FE model Dispersion Solution
% tic
for j1 = 1:n_kap

    tstart = tic;

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
    C_hat = (T_per_k'*Cf*T_per_k);
    M_hat = (T_per_k'*Mf*T_per_k);
    
    %symmetrize matrices
    K_hat = (1/2)*(K_hat+K_hat');
    C_hat = (1/2)*(C_hat+C_hat');
    M_hat = (1/2)*(M_hat+M_hat');
    
    %form symmetric state space matrices
    use_symmetric = false;
%     use_symmetric = true;
    if use_symmetric
            
        % symmetric Equations 2
        A = [zeros(n_dof_per),K_hat;K_hat,C_hat];
        B = [K_hat',zeros(n_dof_per);zeros(n_dof_per),-M_hat];
        A = sparse(A);
        B = sparse(B);        
        A = (1/2)*(A + A');
        B = (1/2)*(B + B');
        

    else

        %form well conditioned state space matrices
        density = 1;
        if density<0.05
            [Krow,Kcol,Kval] = find(K_hat);
            [Mrow,Mcol,Mval] = find(M_hat);
            [Crow,Ccol,Cval] = find(C_hat);

            Arow = [Krow;(1:n_dof_per)'+n_dof_per];
            Acol = [Kcol;(1:n_dof_per)'+n_dof_per];
            Aval = [Kval;ones(n_dof_per,1)];
            A = sparse(Arow,Acol,Aval);

            Brow = [Crow;Mrow;(1:n_dof_per)'+n_dof_per];
            Bcol = [Ccol;Mcol+n_dof_per;(1:n_dof_per)'];
            Bval = [-Cval;-Mval;ones(n_dof_per,1)];
            B = sparse(Brow,Bcol,Bval);

        else
            A = [K_hat,             zeros(n_dof_per);...
                 zeros(n_dof_per),  eye(n_dof_per)];
            B = [-C_hat,            -M_hat;...
                 eye(n_dof_per),    zeros(n_dof_per)];
            A = sparse(A);
            B = sparse(B);
        end
    end

    % Use direct eigenvalue solver if matrices are small enough, otherwise
    % use an iterative solver
    full_eig = 2*n_dof_per<=800  | n_curves>n_dof_per;
    if full_eig
        [PHIs,L_disp] = eig(full(A),full(B),'chol'); 
    else
        [PHIs,L_disp] = eigs(A,B,n_curves,shift);
    end   

    % store eigenvectors and eigenvalues    
    [lambdas,isort] = sort(diag(L_disp));
    lambda(:,j1) = lambdas;
    PHI(:,:,j1) = PHIs(1:n_dof_per,isort);
    
    % display loop timing info
    t_kloop(j1) = toc(tstart);
    fprintf('k-point %i of %i, solution time: %4.2f\n',j1,n_kap,t_kloop(j1))
end

%% sort eigenvalues and eigenvectors
sort_type = 'omega_d';
switch sort_type
    case 'omega_d'
        % natural frequency
        w_n_p = sqrt(real(lambda).^2 + imag(lambda).^2);

        % damping ratio
        zeta = -real(lambda)./w_n_p;

        % damped natural frequency
        omega_d = w_n_p.*sqrt(1-zeta.^2);
%         omega_d = imag(lambda);  


        for i = 1:n_kap
            [omega_d(:,i),isort] = sort(omega_d(:,i));
            lambda(:,i) = lambda(isort,i);
            PHI(:,:,i) = PHI(:,isort,i);
        end
        
        % discard negative frequencies
        disc_neg = false;
        if disc_neg            
            omega_d_avg = sum(omega_d,2)/n_kap;
            tol = max(abs(omega_d(:)))*1e-3;
            i_keep = omega_d_avg>tol;
            lambda = lambda(i_keep,:);
            PHI = PHI(:,i_keep,:);
        end
        
    case 'lambda_abs'
        
        for i = 1:n_kap
            [lambda(:,i),isort] = sort(lambda(:,i));
            PHI(:,:,i) = PHIs(1:n_dof_per,isort);
        end
end