%% Reduced Bloch Mode Expansion dispersion computation Code
% ======================================================================= %
% Dimitri Krattiger
% 1-27-2014

function [omega_RBME,PHI_RBME] = ...
    RBME(K_free,M_free,n_exp,m_exp,KappaPts,i_mode_save,kappa,r1,r2,r3,...
         i_i,i_l,i_r,i_f,i_b,i_d,i_t,i_lf,i_rf,i_rb,i_lb,...
         i_ld,i_rd,i_rt,i_lt,i_fd,i_bd,i_bt,i_ft,...
         i_lfd,i_rfd,i_rbd,i_lbd,i_lft,i_rft,i_rbt,i_lbt)


%% inputs
% K_free = free stiffness matrix
% 
% M_free = free mass matrix
% 
% n_exp = number of expansion points
% 
% m_exp = number of expansion modes
% 
% KappaPts = Create Band Structure by connecting these points in 
%            wave-vector space
%
% i_mode_save = list of modes to keep
%
% kappa
%
% Lx,Ly,Lz,...
%          i_i,i_l,i_r,i_f,i_b,i_d,i_t,i_lf,i_rf,i_rb,i_lb,...
%          i_ld,i_rd,i_rt,i_lt,i_fd,i_bd,i_bt,i_ft,...
%          i_lfd,i_rfd,i_rbd,i_lbd,i_lft,i_rft,i_rbt,i_lbt
     
     
     
% number of curves to save and number of curves to compute
n_curves = numel(i_mode_save);
n_curves_max = max(i_mode_save);    
 
% number of full model DOFs
n_dof = size(K_free,1);

% number of k-points
n_kap = length(kappa);

% % number of high symmetry points
% n_sym_pts = size(KappaPts,2)-1;

% number of Brillouin Zone segments
i_real = ~isnan(KappaPts(1,:));
i_full_seg = find(i_real(1:end-1) & i_real(2:end));
n_segs = sum(i_real(1:end-1) & i_real(2:end));

% % number of dimensions
n_dim = size(kappa,1);

% number of k-points per brillouin zone segment (incl. both end points)
n_kap_seg = 1+(n_kap-1)/n_segs;

% index containing high symmetry points
i_exp = linspace(1,n_kap,(n_exp-1)*n_segs+1);

%% Form periodicity transformation
% ======================================================================= %

% form transformation to enforce bloch boundary conditions
if n_dim == 3
    [T_per0,T_per1,T_per2,T_per3,T_per12,T_per13,T_per23,T_per123]...
        = Periodic_Boundary_Conditions(i_i,i_l,i_r,...
                i_f,i_b,i_d,i_t,...
                i_lf,i_rf,i_rb,i_lb,...
                i_ld,i_rd,i_rt,i_lt,...
                i_fd,i_bd,i_bt,i_ft,...
                i_lfd,i_rfd,i_rbd,i_lbd,...
                i_lft,i_rft,i_rbt,i_lbt);
elseif n_dim == 2
    [T_per0,T_per1,~,T_per2,~,T_per12,~,~]...
        = Periodic_Boundary_Conditions(i_i,...
                         i_l,i_r,...
                         [],[],...
                         i_d,i_t,...
                         [],[],[],[],...
                         i_ld,i_rd,i_rt,i_lt,...
                         [],[],[],[],...
                         [],[],[],[],[],[],[],[]);
end

% periodic model size
n_dof_per = size(T_per0,2);

%% loop through segments of brillouin zone and compute RBME transformation
% ======================================================================= %

% preallocate arrays
omega_RBME = zeros(n_curves,n_kap);
PHI_RBME = zeros(n_dof,n_curves,n_kap);
T_RBME = zeros(n_dof_per,n_exp*m_exp,n_segs);

%initialize counters and  timers
q_RBME = [zeros(1,n_segs*(n_exp-1)),n_kap];
tic
count = 1;
% for i = 1:n_segs
%     for j = 1:n_exp
%         
%         tstart = toc;
% 
%         % index for current k-point
%         q = i_exp((i-1)*(n_exp-1)+j);
%         q_RBME(count) = q;
%         count = count + 1;        
%         
%         if i>1 && j == 1
%             T_RBME(:,1:m_exp,i) = ...
%                 T_RBME(:,(n_exp-1)*m_exp+1:(n_exp)*m_exp,i-1);
%         elseif i == n_segs && j == n_exp && all(KappaPts(:,1) == KappaPts(:,end))
%                 T_RBME(:,(n_exp-1)*m_exp+1:(n_exp)*m_exp,n_segs) = ...
%                     T_RBME(:,1:m_exp,1);
% 
%                 % save full solutions in RBME slots
%                 omega_RBME(:,n_kap) = omega_RBME(:,1);
%                 PHI_RBME(:,:,n_kap) = PHI_RBME(:,:,1);
%                 
%                 disp('added to end')
%         else
% 
% %         lam1  = exp(-1i*kvec'*r1);
% %         lam2  = exp(-1i*kvec'*r2);
% %         if n_dim == 3
% %             lam3  = exp(-1i*kvec'*r3);
% %         end    
%             kvec = kappa(:,q);
% 
%             % assemble periodicity transformation for current k-point
%             if n_dim == 1
%                 % current phase differences
%                 lam1  = exp(-1i*kvec'*r1);
% 
%                 % periodicity transformation
%                 T_per_k = T_per0 + T_per1*lam1;
%             elseif n_dim == 2
% 
%                 % current phase differences
%                 lam1  = exp(-1i*kvec'*r1);
%                 lam2  = exp(-1i*kvec'*r2);
% 
%                 % periodicity transformation
%                 T_per_k = T_per0 + T_per1*lam1 + T_per2*lam2 + T_per12*lam1*lam2;
%             elseif n_dim == 3
% 
%                 % current phase differences
%                 lam1  = exp(-1i*kvec'*r1);
%                 lam2  = exp(-1i*kvec'*r2);
%                 lam3  = exp(-1i*kvec'*r3);
% 
%                 % periodicity transformation
%                 T_per_k = T_per0 + T_per1*lam1 + T_per2*lam2 + T_per3*lam3...
%                     + T_per12*lam1*lam2 + T_per23*lam2*lam3 + T_per13*lam1*lam3...
%                     + T_per123*lam1*lam2*lam3;
%             end           
% 
%             % apply periodicity transformation to mass and stiffness matrices
%             K_hat = (T_per_k'*K_free*T_per_k);
%             M_hat = (T_per_k'*M_free*T_per_k);
% 
%             % apply periodicity transformation to mass and stiffness matrices
%             K_hat = (1/2)*(K_hat'+K_hat);
%             M_hat = (1/2)*(M_hat'+M_hat);    
% 
%             % eigensolution
%             [PHIs,L_disp] = eigs(K_hat,M_hat,m_exp,'sm');
% 
%             % frequency sorting
%             [w_disp,i_disp] = sort((sqrt(diag(L_disp))));
%             PHIs = PHIs(:,i_disp);
%             T_RBME(:,(j-1)*m_exp+1:(j)*m_exp,i) = PHIs;
% 
%             % save full solutions in RBME slots
%             omega_RBME(:,q) = w_disp(i_mode_save);
%             PHI_RBME(:,:,q) = T_per_k*PHIs(:,i_mode_save);
% 
%             disp([n_exp,m_exp,i,j,toc-tstart])
%             
%         end
%     end
% end

% loop through expansion points in each segment
for i = 1:n_segs
    for j = 1:n_exp
        
        tstart = toc;

        % index for current k-point
        q = i_exp((i-1)*(n_exp-1)+j);
        q_RBME(count) = q;
        count = count + 1;        
        
        if i>1 && j == 1 && (i_full_seg(i-1)+1) == i_full_seg(i)
            
            T_RBME(:,1:m_exp,i) = ...
                T_RBME(:,(n_exp-1)*m_exp+1:(n_exp)*m_exp,i-1);
            
        elseif i == n_segs && j == n_exp && all(KappaPts(:,1) == KappaPts(:,end))
            
                T_RBME(:,(n_exp-1)*m_exp+1:(n_exp)*m_exp,n_segs) = ...
                    T_RBME(:,1:m_exp,1);

                % save full solutions in RBME slots
                omega_RBME(:,n_kap) = omega_RBME(:,1);
                PHI_RBME(:,:,n_kap) = PHI_RBME(:,:,1);
                
                disp('added to end')
        else
            
            alpha = (j-1)/(n_exp-1);
            kvec = KappaPts(:,i_full_seg(i))*(1-alpha) + ...
                   KappaPts(:,i_full_seg(i)+1)*alpha;
                        
            % assemble periodicity transformation for current k-point
            if n_dim == 1
                % current phase differences
                lam1  = exp(-1i*kvec'*r1);

                % periodicity transformation
                T_per_k = T_per0 + T_per1*lam1;
            elseif n_dim == 2

                % current phase differences
                lam1  = exp(-1i*kvec'*r1);
                lam2  = exp(-1i*kvec'*r2);

                % periodicity transformation
                T_per_k = T_per0 + T_per1*lam1 + T_per2*lam2 + T_per12*lam1*lam2;
            elseif n_dim == 3

                % current phase differences
                lam1  = exp(-1i*kvec'*r1);
                lam2  = exp(-1i*kvec'*r2);
                lam3  = exp(-1i*kvec'*r3);

                % periodicity transformation
                T_per_k = T_per0 + T_per1*lam1 + T_per2*lam2 + T_per3*lam3...
                    + T_per12*lam1*lam2 + T_per23*lam2*lam3 + T_per13*lam1*lam3...
                    + T_per123*lam1*lam2*lam3;
            end           

            % apply periodicity transformation to mass and stiffness matrices
            K_hat = (T_per_k'*K_free*T_per_k);
            M_hat = (T_per_k'*M_free*T_per_k);

            % apply periodicity transformation to mass and stiffness matrices
            K_hat = (1/2)*(K_hat'+K_hat);
            M_hat = (1/2)*(M_hat'+M_hat);    

            % eigensolution
            [PHIs,L_disp] = eigs(K_hat,M_hat,m_exp,'sm');

            % frequency sorting
            [w_disp,i_disp] = sort((sqrt(diag(L_disp))));
            PHIs = PHIs(:,i_disp);
            T_RBME(:,(j-1)*m_exp+1:(j)*m_exp,i) = PHIs;

            % save full solutions in RBME slots
            omega_RBME(:,q) = w_disp(i_mode_save);
            PHI_RBME(:,:,q) = T_per_k*PHIs(:,i_mode_save);

            disp([n_exp,m_exp,i,j,toc-tstart])
            
        end
    end
end


%% loop through segments of brillouin zone and orthonormalize RBME transformation
for i = 1:n_segs
    [T_RBME(:,:,i),~] = qr(T_RBME(:,:,i),0);
end

%% loop through segments of brillouin zone and compute dispersion
% ======================================================================= %
tic
for i = 1:n_segs
    for j = 1:n_kap_seg
        tstart = toc;

        % k-point index
        q = (i-1)*(n_kap_seg-1) + j;

        if ~ismember(q,q_RBME);
            kvec = kappa(:,q);
            
            % assemble periodicity transformation for current k-point
            if n_dim == 1
                % current phase differences
                lam1  = exp(-1i*kvec'*r1);

                % periodicity transformation
                T_per_k = T_per0 + T_per1*lam1;
            elseif n_dim == 2
                % current phase differences
                lam1  = exp(-1i*kvec'*r1);
                lam2  = exp(-1i*kvec'*r2);

                % periodicity transformation
                T_per_k = T_per0 + T_per1*lam1 + T_per2*lam2 + T_per12*lam1*lam2;
            elseif n_dim == 3
                % current phase differences
                lam1  = exp(-1i*kvec'*r1);
                lam2  = exp(-1i*kvec'*r2);
                lam3  = exp(-1i*kvec'*r3);

                % periodicity transformation
                T_per_k = T_per0 + T_per1*lam1 + T_per2*lam2 + T_per3*lam3...
                    + T_per12*lam1*lam2 + T_per23*lam2*lam3 + T_per13*lam1*lam3...
                    + T_per123*lam1*lam2*lam3;
            end        
            
            % combine periodicity and RBME transformations
            T_per_k_T_RBME = T_per_k*T_RBME(:,:,i);

            % apply periodicity transformation to mass and stiffness matrices
            K_hat = (T_per_k_T_RBME'*K_free*T_per_k_T_RBME);
            M_hat = (T_per_k_T_RBME'*M_free*T_per_k_T_RBME);
            
            % remove zero stiffness nodes
            i_zero = all(K_hat == 0);
            K_hat = K_hat(~i_zero,~i_zero);
            M_hat = M_hat(~i_zero,~i_zero);
            
            % apply periodicity transformation to mass and stiffness matrices
            K_hat = (1/2)*(K_hat'+K_hat);
            M_hat = (1/2)*(M_hat'+M_hat);   
            
            % eigensolution
            [PHIs,L_disp] = eigs(K_hat,M_hat,n_curves_max,'sm');

            % frequency sorting
            [w_disp,i_disp] = sort((sqrt(diag(L_disp))));
            PHIs = PHIs(:,i_disp);
            PHI_RBME(:,:,q) = T_per_k_T_RBME*PHIs(:,i_mode_save);
            omega_RBME(:,q) = w_disp(i_mode_save)';

            % display loop info
            disp([q,toc-tstart])

        end
    end
end