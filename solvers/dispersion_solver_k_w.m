function [kappa,PHI,varargout] = dispersion_solver_k_w(omegas,Kf,Cf,Mf,dof_sets,R,varargin)

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

%% Default options
% ======================================================================= %

% if the last argument given is not an options structure, assume that it 
% is a single parameter specifying the number of desired curves. (This is
% how older versions of the code were called.)
if nargin<6
    options.dummy = 0;
else
    if isstruct(varargin{1})
        options = varargin{1};
    else
        options.n_curves = varargin{1};
    end
end

% does the model have residual enhancement?
resPlus = isstruct(Kf);

% Default behavior for dynamic reduction: 
% If eigenvalue solution method is direct, we really need dynamic
% reduction unless the model is very small. For iterative solution, we are
% taking advantage of sparsity so dynamic reduction is not necessarily
% beneficial.
dynRed = true;
if isfield(options,'full_eig')
    if ~options.full_eig    
        dynRed = false;
    end   
end

% Don't use dynamic reduction if residual-enhanced model is detected. It
% doesn't seem to work. Hopefully model reduction is sufficient so this 
% won't play a large role anyway.
if resPlus    
    dynRed = false;
end

% default option values
defaults.n_curves           = 10;
defaults.verbose            = true;
defaults.full_eig           = true;
defaults.dynamicReduction   = dynRed;
defaults.useConditionValue  = true;

% copy default values into blank fields of options structure
options = setstructfields(defaults,options);

% guess values to be used as centering frequencies by iterative solver
guess1 = 1;
guess2 = 0.0;
guess = guess1;

% Notes about options parameters:
% ===============================
% (1)   Dynamic reduction makes use of a system solution where the matrices
%       are sparse. MATLAB automatically chooses the cholmod when the system
%       to be solved is symmetric positive definite (SPD). When the
%       frequency is zero, the dynamic matrix is SPD, but as the frequency
%       is increased, eventually it becomes indefinite. When this happens
%       a slower solver must be used. Run "spparms('spumoni',1)" before
%       script to output detailed info about the sparse routines being
%       called by MATLAB.
% (2)   The default behavior is to perform full (direct) eigenvalue
%       solutions rather than utilizing an iterative solver. The reason for
%       this is that the eigenvalues are of the form:
%       lambda=exp(i*kappa*Lx). If we find the k smallest lambdas, it does
%       not mean that we have found the k smallest kappas. To guarantee
%       that we have found the desired solutions we need to compute all of
%       the eigenvalues. There may be a way to transform the problem to
%       guarantee this but this is the subject of further research. Note,
%       this would not be a problem if we used the Bloch-operator approach
%       rather than the Bloch boundary condition approach.

%% Prepare for dynamic reduction
% ======================================================================= %
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

n_dof = size(Kf,1);
if isempty(Cf)
    clear Cf;
    if resPlus
        Cf.w0=spalloc(n_dof,n_dof,0);
        Cf.w2=spalloc(n_dof,n_dof,0);
        Cf.w4=spalloc(n_dof,n_dof,0);
    else
        Cf = spalloc(n_dof,n_dof,0);
    end
end

%% Form periodicity transformation
% ======================================================================= %
% form Bloch-periodicity matrices using boundary DOF sets
% Output T_per is a structure with fields s0, s1, s2, s3, s12, s13, s23,
% s123. These can be combined to form the Bloch-periodicity transformation
% matrix for any wavevector. At this point, the solver only uses s0 and s1
% because we only solve in the Gamma-X direction.
if options.dynamicReduction
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
PHI = zeros(n_dof_per,options.n_curves,n_om);
kappa = nan(options.n_curves,n_om);
t_wloop = zeros(1,n_w);
Lsave = nan(options.n_curves,n_om);

% GAMMA-X Direction
Tx = T_per.s1 + T_per.s12 + T_per.s13 + T_per.s123;
T0 = T_per.s0 + T_per.s2 + T_per.s3 + T_per.s23;

% Loop through frequencies to calculate Dispersion Solution
for j1 = 1:n_w

    tstart = tic;
    
    % phase modification at boundaries due to wavevector
    w = omegas(j1);
    
    if resPlus
        
        % residual-enhanced representation
        Kfe = Kf.w0 + w^2*(Kf.w2+(Kf.w2)') + w^4*Kf.w4;
        Cfe = Cf.w0 + w^2*(Cf.w2+(Cf.w2)') + w^4*Cf.w4;
        Mfe = Mf.w0 + w^2*(Mf.w2+(Mf.w2)') + w^4*Mf.w4;
        
    else
        % standard representation
        Kfe = Kf;
        Cfe = Cf;
        Mfe = Mf;
    end
    
    % Form Dynamic matrix
    D = Kfe + 1i*w*Cfe - w^2*Mfe;
    D = (1/2)*(D+D');
    
    % dynamic reduction
    if options.dynamicReduction
        Dhat = D(i_b,i_b)-D(i_b,i_i)*(D(i_i,i_i)\D(i_i,i_b));
        dynRedString = 'dynamically reduced, ';
    else
        Dhat = D;
        dynRedString = [];
    end
        
    % Apply periodicity transformation and decompose D by powers of 
    % e^(i kx Lx)
    D0 = Tx'*Dhat*T0;               % constant
    D1 = T0'*Dhat*T0 + Tx'*Dhat*Tx; % linear in exp(i kx lx)
    D2 = T0'*Dhat*Tx;               % quadratic in exp(i kx lx)
    
    % scaling value to improve conditioning of state space problem. I
    % haven't observe this to make a noticeable difference.
    if options.useConditionValue
        condVal = (sum(diag(D0))+sum(diag(D1))+sum(diag(D2)))/(3*n_dof_per);
    else
        condVal = 1;
    end
        
    % prepare for sparse represntation by finding indices of non-zero
    % entries
    [D0row,D0col,D0val] = find(D0);
    [D1row,D1col,D1val] = find(D1);
    [D2row,D2col,D2val] = find(D2);

    % sparse representation of state space matrices
    Arow = [D0row;((n_dof_per+1):(2*n_dof_per))'];
    Acol = [D0col;((n_dof_per+1):(2*n_dof_per))'];
    Aval = [D0val;ones(n_dof_per,1)*condVal];

    Brow = [D1row;D2row;((n_dof_per+1):(2*n_dof_per))'];
    Bcol = [D1col;D2col+n_dof_per;(1:n_dof_per)'];
    Bval = [-D1val;-D2val;ones(n_dof_per,1)*condVal];

    % form sparse state space matrices
    A = sparse(Arow,Acol,Aval,(2*n_dof_per),(2*n_dof_per));
    B = sparse(Brow,Bcol,Bval,(2*n_dof_per),(2*n_dof_per));

    % % form state space matrices
    % A = [D0,zeros(n_dof_per);zeros(n_dof_per),eye(n_dof_per)*condVal];
    % B = [-D1,-D2;eye(n_dof_per)*condVal,zeros(n_dof_per)];
    
    % eigenvalue solution
    if options.full_eig
        [PHIs,L] = eig(full(A),full(B),'vector');
        
        % solution string
        soln_type = 'direct';
    else
        
        % The matrices A and B should be sparse
%         try
            [PHIs,L] = eigs((A),(B),options.n_curves,guess);
            L = diag(L);
            
%         catch
%             
%             % note that inverse arnoldi can fail if we are looking for
%             % solutions near a pt that happens to exactly contain a
%             % solution. Thus, if the iterative solver fails, we try
%             % perturbing the centering point.
%             if guess == guess1
%                 guess = guess2;
%                 oldguess = guess1;
%             else
%                 guess = guess1;
%                 oldguess = guess2;
%             end
%             
%             fprintf(['"eigs" solver failed looking for solutions near ',...
%                     num2str(oldguess),'. Trying to find solutions near ',...
%                     num2str(guess),'.\n'])
%                 
%             [PHIs,L] = eigs(A,B,options.n_curves,guess);
%             L = diag(L);
%         end
        
        % solution string (used for outputing info to command window)
        soln_type = 'iterative';
    end
    
    % post process solutions     
    % gamma-X direction
    logL = log(L');
    kappas = logL/(-1i*R(1));
    
    % sort solutions by complex magnitude
    [~,i_sort] = sort(abs(kappas));
    kappas = kappas(i_sort);
    PHIs = PHIs(1:n_dof_per,i_sort);
    L = L(i_sort);
    
    % Save solutions
    PHI(:,:,j1) = PHIs(:,1:options.n_curves);
    kappa(:,j1) = kappas(1:options.n_curves)';
    Lsave(:,j1) = L(1:options.n_curves);
    
    % loop timing
    t_wloop(j1) = toc(tstart);
    
    % display loop timing info
    fprintf(['freq. point %i of %i, ',...
             soln_type,' solution, ',...
             dynRedString,...
             'calc. time: %4.2f\n'],j1,n_w,t_wloop(j1))
end

% if enough arguments are requested, output the w-loop timing vector
if nargout>=3
    varargout{1} = t_wloop;
end