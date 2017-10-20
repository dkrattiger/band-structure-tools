%% Group Velocity Direct Computation Code
% Dimitri Krattiger
% 4-9-2014

function v_g = group_velocity(K_free,M_free,PHI,omega,n_curves,kappa,Lx,Ly,Lz,...
         dof_sets)

     
 t_start_vg = tic;
     
% number of dimensions
n_dim = size(kappa,1);

% number of full model DOFs
n_dof = size(K_free,1);

% number of k-points
n_kap = length(kappa);

%% Form Periodicity Transformation Matrix
% ======================================================================= %
[T_per0,T_perx,T_pery,T_perz,T_perxy,T_perxz,T_peryz,T_perxyz]...
    = Periodic_Boundary_Conditions_sparse(dof_sets);

%% compute group velocities
% ======================================================================= %

% preallocate arrays
v_g = zeros(n_curves,n_kap);
tol = (1e-5)*max(max(omega));

for j2 = 1:n_kap
    j2
    t_kloop_start = toc(t_start_vg);

    % wave vector at current k-point
    kx = kappa(1,j2);
    ky = kappa(2,j2);
    if n_dim == 3
        kz = kappa(3,j2);
    end

    % phase modification at boundaries due to wavevector
    lamx = exp(-1i*kx*Lx);
    lamy = exp(-1i*ky*Ly);
    if n_dim == 3
        lamz = exp(-1i*kz*Lz);
    end
    
    % direction of differential k-vector
    if j2>1
        k_hat = kappa(:,j2)-kappa(:,j2-1);
        k_hat = k_hat/norm(k_hat);
    else
        k_hat = kappa(:,j2+1)-kappa(:,j2);
        k_hat = k_hat/norm(k_hat);
    end
    
    % form bloch boundary condition transformation at current k-point
    if n_dim == 1
        T_per_k = T_per0 + T_perx*lamx;
    elseif n_dim == 2
        T_per_k = T_per0 + T_perx*lamx + T_pery*lamy + T_perxy*lamx*lamy;
    elseif n_dim == 3
        T_per_k = T_per0 + T_perx*lamx + T_pery*lamy + T_perz*lamz...
            + T_perxy*lamx*lamy + T_peryz*lamy*lamz + T_perxz*lamx*lamz...
            + T_perxyz*lamx*lamy*lamz;
    end
    
    % form mass matrix
    M = T_per_k'*M_free*T_per_k;  
    M = (1/2)*(M+M');
%     K = T_per_k'*K_free*T_per_k; 
%     K = (1/2)*(K+K');
    
    
%     [PHI_k,L] = eig(full(K),full(M));
%     [omega_k,i_sort] = sort(sqrt(diag(L)));
%     PHI_k = PHI_k(:,i_sort);
    
    % mass normalize modes
    PHI_k = PHI(:,:,j2);
    omega_k = omega(:,j2);
    for i = 1:n_curves
        PHI_k(:,i) = PHI_k(:,i)/sqrt(PHI_k(:,i)'*M*PHI_k(:,i));
    end
    
%     PHI_k(:,[1:10])'*M*PHI_k(:,[1:10])
%     pause
    
    % form bloch boundary condition transformation at current k-point
    if n_dim == 1
        
    elseif n_dim == 2
        
        % periodicity transformation derivative
        dTdkx = -1i*Lx*T_perx*lamx - 1i*Lx*T_perxy*lamx*lamy; 
        dTdky = -1i*Ly*T_pery*lamy - 1i*Ly*T_perxy*lamx*lamy;         
        
        K_prime_x = dTdkx'*K_free*T_per_k + T_per_k'*K_free*dTdkx;
        K_prime_y = dTdky'*K_free*T_per_k + T_per_k'*K_free*dTdky;
        K_prime = K_prime_x*k_hat(1) + K_prime_y*k_hat(2);
        
        M_prime_x = dTdkx'*M_free*T_per_k + T_per_k'*M_free*dTdkx;
        M_prime_y = dTdky'*M_free*T_per_k + T_per_k'*M_free*dTdky;
        M_prime = M_prime_x*k_hat(1) + M_prime_y*k_hat(2);
        
    elseif n_dim == 3
        % periodicity transformation derivative
        dTdkx = -1i*Lx*(T_perx*lamx + T_perxy*lamx*lamy + T_perxz*lamx*lamz + T_perxyz*lamx*lamy*lamz);
        dTdky = -1i*Ly*(T_pery*lamy + T_perxy*lamx*lamy + T_peryz*lamy*lamz + T_perxyz*lamx*lamy*lamz);
        dTdkz = -1i*Lz*(T_perz*lamz + T_perxz*lamx*lamz + T_peryz*lamy*lamz + T_perxyz*lamx*lamy*lamz); 
        
        K_prime_x = dTdkx'*K_free*T_per_k + T_per_k'*K_free*dTdkx;
        K_prime_y = dTdky'*K_free*T_per_k + T_per_k'*K_free*dTdky;
        K_prime_z = dTdkz'*K_free*T_per_k + T_per_k'*K_free*dTdkz;
        K_prime = K_prime_x*k_hat(1) + K_prime_y*k_hat(2) + K_prime_z*k_hat(3);
        
        M_prime_x = dTdkx'*M_free*T_per_k + T_per_k'*M_free*dTdkx;
        M_prime_y = dTdky'*M_free*T_per_k + T_per_k'*M_free*dTdky;
        M_prime_z = dTdkz'*M_free*T_per_k + T_per_k'*M_free*dTdkz;
        M_prime = M_prime_x*k_hat(1) + M_prime_y*k_hat(2) + M_prime_z*k_hat(3);
    end
    
    K_prime = (1/2)*(K_prime+K_prime');
    M_prime = (1/2)*(M_prime+M_prime');
    
    % compute group velocity    
    i = 1;
    while i<=n_curves
        
        % index of eigenvalues that are the same as omega_i
        i_curves = find(abs(omega_k-omega_k(i))<=tol);
        n_deg = length(i_curves);
        
        % ith eigenvalue
        w = omega_k(i);  
%         [omega(i_curves),omega(i_curves).^2];
        
        % degenerate mode shapes
        phi = PHI_k(:,i_curves);
        
        % group velocity
        if n_deg == 1
%             v_g(i,j2) = (1/(2*w))*phi'*(K_prime-w^2*M_prime)*phi;
            v_g(i,j2) = phi'*(K_prime-w^2*M_prime)*phi;
        else
%             phi
            A = phi'*(K_prime-w^2*M_prime)*phi;
            A = (1/2)*(A+A');
%             v_g(i_curves,j2) = (1/(2*w))*(eig(A));
            v_g(i_curves,j2) = (eig(A));
%             w
%             eig(A)
%             pause
        end
        
%         pause
        i = i + n_deg;
    end
end

