function [omega,PHI,t_kloop,varargout] = dispersion_solver_w_k(kappa,K,M,dof_sets,R,varargin)

% Dimitri Krattiger
%
% Description
% ===========
% This code computes the band-structure frequencies at each wave-vector
% specified in "kappa".
%
% inputs
% ======
% kappa     = list of wave vectors at which to compute dispersion
%             frequencies
%
% K         = Free stiffness matrix. Can be a reduced representation (e.g.
%             BMS)
% 
% M         = Free mass matrix. Can be a reduced representation (e.g. BMS)
%
% dof_sets  = structure containing indices of node sets (obtain with 
%             "find_node_sets.m" and "node2dof.m" functions)
%
% R         = Lattice vectors: R = [r1,r2,r3]
%
% n_curves  = number of band-structure curves to compute
%
% outputs
% =======
% omega     	= band-structure frequencies
% 
% PHI       	= band-structure mode shapes
%
% t_kloop       = computation time for each k point


%% Add subfolders with dependent libraries to Matlab path
% ======================================================================= %
addpath(genpath('libraries'))

%% Default options
% ======================================================================= %

defaults.n_curves = 10;
defaults.wCenter  = 0;
defaults.verbose  = true;
defaults.k        = 30;
% defaults.storephi = n_dof_per<10000 ;

%defaults.fullEig = gets determined by model size and # curves;
%defaults.fullIterativeEig = gets determined by model size and # curves;

if nargin<6
    options.dummy = 0;
else
    if isstruct(varargin{1})
        options = varargin{1};
    else
        options.n_curves = varargin{1};
    end
end
options = setstructfields(defaults,options);

% verbose = true;
resPlus = isstruct(K);


%% Sparsify matrices if density is less than 0.2
% ======================================================================= %

% value of matrix density below which to use sparse algorithms
density_tol = 0.1; % 0.1 is quite high, not sure what the exact cutoff should be though 
if resPlus
    density = nnz(M.w0)/numel(M.w0);
    if density < density_tol
        K.w0 = sparse(K.w0);        M.w0 = sparse(M.w0);
        K.w2 = sparse(K.w2);        M.w2 = sparse(M.w2);
        K.w4 = sparse(K.w4);        M.w4 = sparse(M.w4);
    else
        K.w0 = full(K.w0);        M.w0 = full(M.w0);
        K.w2 = full(K.w2);        M.w2 = full(M.w2);
        K.w4 = full(K.w4);        M.w4 = full(M.w4);
    end    
else
    density = nnz(K)/numel(K);
    if density < density_tol
        K = sparse(K);        M = sparse(M);
    else
        K = full(K);          M = full(M);
    end
end
[~,n_kap] = size(kappa);

%% Form periodicity transformation
% ======================================================================= %
% form Bloch-periodicity matrices using boundary DOF sets
% Output T_per is a structure with fields s0, s1, s2, s3, s12, s13, s23,
% s123. These can be combined to form the Bloch-periodicity transformation
% matrix for any wavevector.
T_per = Periodic_Boundary_Conditions(dof_sets);

% number of DOFs
n_dof = size(T_per.s0,1);
n_dof_per = size(T_per.s0,2);


%% Update computed default options if not given
% ======================================================================= %
% number of dispersion curves to compute
n_curves = options.n_curves;
      
if options.verbose
    fprintf('Model Size: %iDOF (%i before enforcing periodicity)\n',...
    n_dof_per,n_dof);
end

% update full eigenvalue solution options
if ~isfield(options,'fullEig')
    options.fullEig = n_dof_per/n_curves<=5  || n_dof_per<1000; %n_dof_per<800;
end
if ~isfield(options,'fullIterativeEig')
    options.fullIterativeEig = n_dof_per>5000; % rough cutoff.
end

% check if mode shapes should be stored and update that option field
if ~isfield(options,'storePHI')
    options.storePHI = n_dof_per<5000 | ~options.fullEig;
end


%% Dispersion Solution
% ======================================================================= %

% preallocate Arrays
omega = zeros(n_curves,n_kap);
iWindow = zeros(n_curves,n_kap);

% Preallocate mode shape matrix
if options.storePHI
    PHI = zeros(n_dof_per,n_curves,n_kap);
else
    PHI = [];
end    
t_kloop = zeros(1,n_kap);

% Full FE model Dispersion Solution
for j1 = 1:n_kap

    tstart = tic;

    % wave vector at current k-point
    kvec = kappa(:,j1);
    
    % phase modification at boundaries due to wavevector
    lam  = exp(-1i*kvec'*R);
    lam(end+1:3) = 0;
    
    T_per_k = T_per.s0 + T_per.s1*lam(1) + T_per.s2*lam(2) + T_per.s3*lam(3)...
            + T_per.s12*lam(1)*lam(2) + T_per.s23*lam(2)*lam(3) + T_per.s13*lam(1)*lam(3)...
            + T_per.s123*lam(1)*lam(2)*lam(3);
                
    % different approach if residual enhancement
    if resPlus
        
        % Apply Periodicity Transformation to all matrices
        Kp.w0 = (T_per_k'*K.w0*T_per_k);
        Kp.w2 = (T_per_k'*K.w2*T_per_k);
        Kp.w4 = (T_per_k'*K.w4*T_per_k);

        Mp.w0 = (T_per_k'*M.w0*T_per_k);
        Mp.w2 = (T_per_k'*M.w2*T_per_k);
        Mp.w4 = (T_per_k'*M.w4*T_per_k);
        
        % to be used in approximating un-enhanced frequencies
        MinvK = Mp.w0\Kp.w0;

        reduced = true;
        if reduced
            % these matrices are somewhat easier to compute and can be
            % shown to be equivalent to the straightforward versions. More
            % details in:
            %
            % "A simplified error estimator for the CB method and its 
            % application to error control"
            % S. Boo, J.G. Kim, P.S. Lee
            M_hats.w0 = Mp.w0;
            M_hats.w2 = Mp.w2'*MinvK;
            M_hats.w4 = MinvK'*Mp.w4*MinvK;

            K_hats.w0 = Kp.w0;
            K_hats.w4 = MinvK'*M_hats.w2;

            K_hat = K_hats.w0 + K_hats.w4;
            M_hat = M_hats.w0 + M_hats.w2 + M_hats.w2' + M_hats.w4;
            
        else
            
            % straightforward matrices
            M_hats.w0 = Mp.w0;
            M_hats.w2 = Mp.w2'*MinvK;
            M_hats.w4 = MinvK'*Mp.w4*MinvK;

            K_hats.w0 = Kp.w0;
            K_hats.w2 = Kp.w2'*MinvK;
            K_hats.w4 = MinvK'*Kp.w4*MinvK;

            K_hat = K_hats.w0 + K_hats.w2 + K_hats.w2' + K_hats.w4;
            M_hat = M_hats.w0 + M_hats.w2 + M_hats.w2' + M_hats.w4;
        end
    else

        % apply periodicity transformation to mass and stiffness matrices
        K_hat = (T_per_k'*K*T_per_k);
        M_hat = (T_per_k'*M*T_per_k);
        
    end
    
    %symmetrize matrices
    K_hat = (1/2)*(K_hat+K_hat');
    M_hat = (1/2)*(M_hat+M_hat');
        
    % Use direct eigenvalue solver if matrices are small enough, otherwise
    % use an iterative solver
    if options.fullEig
        if options.fullIterativeEig
            [PHIs,L_disp,iWindow(:,j1)] = eigsFull(K_hat,M_hat,options); 
            
        else
            [PHIs,L_disp] = eig(full(K_hat),full(M_hat),'vector'); 
        end
    else
        options_eigs.p = 3*n_curves;
        [PHIs,L_disp] = eigs(K_hat,M_hat,n_curves,options.wCenter^2,options_eigs);
        L_disp = diag(L_disp);
    end

    % frequency sorting
    w_disp = sqrt(L_disp);
    
    % sort dispersion curves by closeness to centering frequency 
    [~,i_disp] = sort(abs(w_disp.^2-options.wCenter^2));

    % truncate number of dispersion frequencies
    w_disp = w_disp(i_disp(1:n_curves));
    PHIs = PHIs(:,i_disp(1:n_curves));
    
    % re-sort dispersion curves by magnitude
    [~,i_disp] = sort(w_disp);
    w_disp = w_disp(i_disp);
    PHIs = PHIs(:,i_disp);
    iWindow(:,j1) = iWindow(i_disp,j1);
    
    % store results
    omega(:,j1) = w_disp;
    if options.storePHI
        PHI(:,:,j1) = PHIs(:,1:n_curves);
    end

    t_kloop(j1) = toc(tstart);
    
    % display loop timing info
    if options.fullEig
        if options.fullIterativeEig
            fprintf('k-point %i of %i, solution time: %4.2f (full iterative solver)\n',...
                j1,n_kap,t_kloop(j1))
        else
            fprintf('k-point %i of %i, solution time: %4.2f (full direct solver)\n',...
                j1,n_kap,t_kloop(j1))
        end
    else
        fprintf('k-point %i of %i, solution time: %4.2f (partial iterative solver)\n',j1,n_kap,t_kloop(j1))
    end

end

%% output 
% ======================================================================= %

if nargout>=4
    varargout{1} = iWindow;
end

%% Remove subfolders from Matlab path
% ======================================================================= %
rmpath(genpath('libraries'))