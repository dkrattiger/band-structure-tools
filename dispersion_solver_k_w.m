function [kappa,PHI,varargout] = dispersion_solver_k_w(omegas,Kf,Cf,Mf,dof_sets,R,n_curves)

% Dimitri Krattiger
%
% DESCRIPTION
% -----------
% This function computes the k(w) dispersion using the
% Bloch-boundary-condition approach. The code currently only supports
% calculation of dispersion in the Gamma-X direction.
%
% INPUTS
% ------
%
% omegas        = vector of frequencies
%
% Kf            = free stiffness matrix (no boundary conditions applied)
%
% Mf            = free mass matrix
%
% Cf            = free damping matrix (use scalar Cf = 0 for no damping)
%
% dof_sets      = structure containing Degree of Freedom sets for each boundary
%                 set. (Use "find_node_sets", "node2dof" functions to
%                 generate this)
%
% R             = Lattice vectors, e.g. R = [[Lx;0;0],[0;Ly;0]] gives
%                 lattice vectors for a 3D model with 2D periodicity
%
% OUTPUTS
% -------
%
% kappa         = array containing wave-vector solutions (n_DOF x n_freq)
%
% PHI           = array containing eigenvector solutions 
%                 (n_DOF x n_DOF x n_freq)
%
% t_wloop       = timing results for each frequency point (n_freq x 1)

full_eig = false;
if isempty(n_curves)
    full_eig = true;
elseif n_curves<0
    full_eig = true;
end

% guess1 = 'sm';
guess1 = 1.1;
guess2 = 0.;
guess = guess1;

% use dynamic reduction?
dyn_red = true;
if ~full_eig    
    dyn_red = false;
end

%% Prepare for dynamic reduction
% ======================================================================= %
% form Bloch-periodicity matrices using boundary DOF sets
% Output T_per is a structure with fields s0, s1, s2, s3, s12, s13, s23,
% s123. These can be combined to form the Bloch-periodicity transformation
% matrix for any wavevector.
i_i = dof_sets.i;

% create a list of all boundary DOFs
dof_set_names = fieldnames(dof_sets);
i_b = [];
for i = 2:length(dof_set_names)
    i_b = [i_b,dof_sets.(dof_set_names{i})'];
end
i_b = sort(i_b);
n_b = length(i_b);

%% Check whether damping matrix is included
% ======================================================================= %
use_C = true;
if isempty(Cf)
    use_C = false;
end

%% Form periodicity transformation
% ======================================================================= %
% form Bloch-periodicity matrices using boundary DOF sets
% Output T_per is a structure with fields s0, s1, s2, s3, s12, s13, s23,
% s123. These can be combined to form the Bloch-periodicity transformation
% matrix for any wavevector.
if dyn_red
    dof_sets.i = [];
end
T_per = Periodic_Boundary_Conditions(dof_sets);

%% Loop through Frequencies
% ======================================================================= %
% number of frequencies, wave propagation dimensions, and 
n_w = length(omegas);
n_p = size(R,2);
[n_dof,n_dof_per] = size(T_per.s0);
n_om = length(omegas);
        
% preallocate arrays
if full_eig
    n_curves = 2*n_dof_per;
end
PHI = zeros(n_dof_per,n_curves,n_om);
kappa = nan(n_curves,n_om);
t_wloop = zeros(1,n_w);
Lsave = nan(n_curves,n_om);

% GAMMA-X Direction
Tx = T_per.s1 + T_per.s12 + T_per.s13 + T_per.s123;
T0 = T_per.s0 + T_per.s2 + T_per.s3 + T_per.s23;

% Loop through frequencies to calculate Dispersion Solution
for j1 = 1:n_w

    tstart = tic;
    
    % phase modification at boundaries due to wavevector
    w = omegas(j1);
    
    % Form Dynamic matrix
    if use_C
        D = Kf+1i*w*Cf-w^2*Mf;
    else
        D = Kf-w^2*Mf;
    end
    
    % dynamic reduction
    if dyn_red
        Dhat = D(i_b,i_b)-D(i_b,i_i)*(D(i_i,i_i)\D(i_i,i_b));
    else
        Dhat = D;
    end
        
    % Apply periodicity transformation and decompose D by powers of 
    % e^(i kx Lx)
    D0 = Tx'*Dhat*T0;
    D1 = T0'*Dhat*T0 + Tx'*Dhat*Tx;
    D2 = T0'*Dhat*Tx;    
    
    [D0row,D0col,D0val] = find(D0);
    [D1row,D1col,D1val] = find(D1);
    [D2row,D2col,D2val] = find(D2);
    
    % sparse representation of state space matrices
    Arow = [D0row;((n_dof_per+1):(2*n_dof_per))'];
    Acol = [D0col;((n_dof_per+1):(2*n_dof_per))'];
    Aval = [D0val;ones(n_dof_per,1)];
    
    Brow = [D1row;D2row;((n_dof_per+1):(2*n_dof_per))'];
    Bcol = [D1col;D2col+n_dof_per;(1:n_dof_per)'];
    Bval = [-D1val;-D2val;ones(n_dof_per,1)];
    
    % form sparse state space matrices
    A = sparse(Arow,Acol,Aval,(2*n_dof_per),(2*n_dof_per));
    B = sparse(Brow,Bcol,Bval,(2*n_dof_per),(2*n_dof_per));
    
    % eigenvalue solution
    if full_eig
        [PHIs,L] = eig(full(A),full(B),'vector');
    else
        try
            [PHIs,L] = eigs((A),(B),n_curves,guess);L = diag(L);
%             disp('iterative, guess = ');disp(guess)
        catch
            if guess == guess1
                guess = guess2;
            else
                guess = guess1;
            end
            [PHIs,L] = eigs((A),(B),n_curves,guess);L = diag(L);
            disp('iterative, guess = ');disp(guess)
        end
%             disp('iterative, guess = ');disp(guess)
%             [PHIs,L] = eigs((A),(B),n_curves,guess);L = diag(L);
            
        
    end
    
% full_eig

    % post process solutions     
    % gamma-X direction
    L_tol = 1e-9;
    L_tol = 1e-99;
%     L(abs(L)<L_tol) = nan + 1i*nan;
    logL = log(L');
    Lsave(:,j1) = L;
    
%     vs = [logL;zeros(n_p-1,length(L))]/(-1i);
%     kapvecs = R(1:n_p,1:n_p)\vs;
%     kappas = kapvecs(1,:);
    kappas = logL/(-1i*R(1));
    
%     if rem(j1,20) == 0
%         figure(212);clf
%         subplot(1,2,1)
%         plot3(real(Lsave),imag(Lsave),omegas,'r.')
%         subplot(1,2,2)
%         plot3(real(kappas),imag(kappas),omegas,'b.')
%         drawnow
%         pause(0.01)
%     end
    % remove infinite kappa values by setting to nan
%     kappas(isinf(abs(kappas))) = nan + 1i*nan;
    
    % sort solutions by complex magnitude
    [~,i_sort] = sort(abs(kappas));
    kappas = kappas(i_sort);
    PHIs = PHIs(1:n_dof_per,i_sort);
%     PHIs = PHIs(1:n_dof_per,1:n_curves);
% size(PHIs)
% size(PHI)
% n_dof_per
% n_dof
    
    % Save solutions
    PHI(:,:,j1) = PHIs;
    kappa(:,j1) = kappas';
    lam(:,j1) = (L);
    
%     if rem(j1,20) == 0
%         figure(212);clf
%         subplot(1,2,1)
%         plot3(real(Lsave),imag(Lsave),omegas,'r.')
%         xlim([-3,3]);ylim([-3,3])
%         subplot(1,2,2)
%         
%         plot3(real(kappa),imag(kappa),omegas,'b.')
% %         xlim([-3,3]);ylim([-3,3])
%         drawnow
%         pause(0.01)
%         pause
%         
%     end   
    
    % loop timing
    t_wloop(j1) = toc(tstart);    
    
    % display loop timing info
    fprintf('freq. point %i of %i, solution time: %4.2f\n',j1,n_w,t_wloop(j1))
end

if nargout>=3
    varargout{1} = t_wloop;
end