function [K_BMS,M_BMS,dof_sets_mod,infoBMS,varargout] = ...
                BMS(K_free,M_free,X,R,options)

            
% Dimitri Krattiger
%
% Sample Call
% ===========
% options.n_FI = 10; % use 10 fixed-interface modes
% options.n_CC = 25; % begin interface reduction with 25 modes
% [K_BMS,M_BMS,dof_sets_mod,t_up_front] = BMS(K_free,M_free,X,R,options);
%
% Description
% ===========
% This code uses the Bloch Mode Synthesis (BMS) model reduction to form a
% reduced order model that can be solved for frequencies at any point in
% the Brillouin zone. The model order reduction is performed on a "free" 
% unit cell model.
%
% inputs
% ======
% K_free    = free stiffness matrix
% 
% M_free    = free mass matrix
%
% X         = nodal coordinates for FEM model 
%            (X=[x(:)] in 1D models
%             X=[x(:),y(:)] in 2D models
%             X=[x(:),y(:),z(:)] in 3D models)
% 
% R         = lattice vectors in 1, 2, and 3 directions 
%             (R = r1 for 1 direction of periodicity i.e. 1P
%              R = [r1,r2] for 2 directions of periodicity i.e. 2P
%              R = [r1,r2,r3] for 3 directions of periodicity i.e. 3P)
%
% options   = structure containing options parameters. Some options must
%             be specified:
%
%                1) options.n_FI -or- options.w_i
%                2) options.n_CC -or- options.w_b  
%                    (unless options.BoundaryMethod = 'none')
%
%             STILL NEED TO ELABORATE ON ADDITIONAL OPTIONS 
%
% outputs
% =======
% K_BMS         = free BMS stiffness matrix
% 
% M_BMS         = free BMS mass matrix
%
% dof_sets_mod  = modified DOF set structure
% 
% t_up_front    = BMS model reduction calculation time (note this does not
%                 include the actual dispersion calculation)
%
% T_BMS         = transformation between full DOF vector and BMS DOF vector

%% Check what inputs are given and set rest to default
% ======================================================================= %

% set_default values for options
defaults.InteriorMethod = 'CB+';
defaults.BoundaryMethod = 'exact';
defaults.n_FI           = [];
defaults.n_CC           = [];
defaults.w_i            = [];
defaults.w_b            = [];
defaults.wQS            = 0;
defaults.verbose        = true;
defaults.verboseTab     = '%%%%  ';
defaults.verboseTabSum  = '';
defaults.plots          = true;
defaults.outputT        = nargout>=5;

if ~exist('options')
    options.dummy = 0;
end

% default to regular CB with no interface reduction
% if a non-zero quasi-static frequency is specified
if isfield(options,'wQS')
    if all(options.wQS~=0)
        defaults.InteriorMethod = 'CB';
        defaults.BoundaryMethod = 'none';
    end
end

% fill in unspecified options fields with the default values
options = setstructfields(defaults,options);

% check if property arguments are acceptable
options = CheckInputs(options);


%% Display Timing info
% ======================================================================= %

t_start_BMS = tic;  
if options.verbose  
    tb = options.verboseTabSum;
    options.verboseTabSum = [options.verboseTabSum,options.verboseTab];
    tbi = options.verboseTabSum;
    fprintf([tb,'\n'])
    fprintf([tb,'Bloch Mode Synthesis (BMS) Reduction\n'])
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
end

%% Code Initial Setup
% ======================================================================= %

% number of full model DOFs
n_dof = size(K_free,1);

% number of nodes
[n_nodes,n_dim] = size(X);

% number of degrees of freedom per node
n_dpn = n_dof/n_nodes;

% sparsify K_free and M_free
K_free = sparse(K_free);
M_free = sparse(M_free);

%% Use geometry to create indices for interior and boundary DOF sets
% ======================================================================= %

% use node coordinates and lattice vectors to get indices of node sets 
% that are on boundaries. Output is a structure containing fields for each
% boundary set
node_sets = find_node_sets(X,R);

% convert node index sets to DOF index sets. Output is a structure
% containing same fields as "node_sets" structure.
dof_sets = node2dof(node_sets,n_dpn);

% take index of all interior DOFs from "dof_sets" structure
i_i_nodes = node_sets.i;

% create a list of all boundary nodes
node_set_names = fieldnames(node_sets);
i_b_nodes = [];
for i = 2:length(node_set_names)
    i_b_nodes = [i_b_nodes,node_sets.(node_set_names{i})];
end

% create DOF lists from node lists (assuming node by node DOF sort)
i_i = node2dof(i_i_nodes,n_dpn);
i_b = node2dof(i_b_nodes,n_dpn);

% length of boundary and interior partitions
n_b = length(i_b);
n_i = length(i_i);

%% Partition Unit Cell if using AMLS (Automated Multi-Level Substructuring)
% ======================================================================= %

if any(strcmpi(options.InteriorMethod,{'amls','amls+'}))
    
    % partition unit cell (boundary DOF sorting will be by set (i.e.
    % l, r, d, t, f, b, etc.)
    max_ss_dofs = 400;          % maximum substructure dofs
    max_ss_nodes = floor(max_ss_dofs/n_dpn);
    [i_ss_dofs] = partition_unit_cell(X,R,K_free,...
        i_i_nodes,i_b_nodes,max_ss_nodes,options.plots);
    
end

%% use AMLS to reduce model, otherwise use CB
% ======================================================================= %

if any(strcmpi(options.InteriorMethod,{'amls','cb'}))
    switch lower(options.InteriorMethod)    
        case 'amls'
            
            % use automated multi-level substructuring to reduce free matrices   
            [K_BMS,M_BMS,dof_sets_mod,T_BMS] = AMLS(i_ss_dofs,K_free,M_free,options.w_i,dof_sets,options);
            
        case 'cb'
            
            % perform CB reduction to reduce free matrices
            [K_BMS,M_BMS,dof_sets_mod,T_BMS,~,~,~,~,wQS] = CB(K_free,M_free,dof_sets,options);     
            
    end

    %% Perform LIR reduction
    t_start_LIR = tic;
    
    if ~strcmpi(options.BoundaryMethod,'none')
        [K_BMS,M_BMS,dof_sets_mod,T_LIR] = LIR(K_BMS,M_BMS,dof_sets_mod,options);
        T_BMS = T_BMS*T_LIR;
    end
    
elseif any(strcmpi(options.InteriorMethod,{'amls+','cb+'}))
    
    switch lower(options.InteriorMethod)
        
        % use residual-enhanced automated multi-level substructuring to reduce free matrices  
        case 'amls+'
            [K_BMS,M_BMS,dof_sets_mod,T_BMS] = ...
                AMLS_plus_LIR(i_ss_dofs,dof_sets,K_free,M_free,options.w_i,options);
            
        % use residual-enhanced HCB to reduce free matrices  
        case 'cb+'
            [K_BMS,M_BMS,dof_sets_mod,T_BMS] = ...
                CB_plus_LIR(K_free,M_free,dof_sets,options);
    end    
end


%% output transformation matrix
% ======================================================================= %
if options.outputT
    varargout{1} = T_BMS;
end

if ~isequal(options.wQS,0) && nargout >=6
    varargout{2} = wQS;
end
    

%% Timing results
% ======================================================================= %
t_up_front = toc(t_start_BMS);

if options.verbose  

    n_FI = length(dof_sets_mod.i);
    if isstruct(K_BMS)
        n_b2 = size(K_BMS.w0,1)-n_FI;
    else
        n_b2 = size(K_BMS,1)-n_FI;
    end

    % display up front timing info
    fprintf([tbi,'No. interior DOF reduced from %i to %i,\n'],n_i,n_FI)
    fprintf([tbi,'No. boundary DOF reduced from %i to %i,\n'],n_b,n_b2)
    fprintf([tbi,'BMS model computation time: %4.2f\n'],t_up_front)
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])  
end

%% Assign variables to info structure for output
% ======================================================================= %

infoBMS.t_up_front = t_up_front;
infoBMS.n_FI = n_FI;
infoBMS.n_b = n_b2;


function opts = CheckInputs(opts)

%% Check if Methods are valid
% ======================================================================= %
% check if interior method is valid
if ~any(strcmpi(opts.InteriorMethod,{'CB','CB+','AMLS','AMLS+'}))
    error('InteriorMethod must be "CB", "CB+", "AMLS" or "AMLS+"');
end

% check if boundary method is valid
if ~any(strcmpi(opts.BoundaryMethod,{'exact','weak','hybrid','none'}))
    error('BoundaryMethod must be "exact", "weak", or "none"');
end

%% perform checks for Interior reduction parameters
% ======================================================================= %

% Cannot do AMLS without an interior cutoff frequency
if any(strcmpi(opts.InteriorMethod,{'amls','amls+'})) && isempty(opts.w_i)
    error('AMLS methods requires an interior cutoff frequency, "w_i", to be specified');
end

% neither w_i nor n_FI is given
if isempty(opts.w_i) && isempty(opts.n_FI)
    error('Either w_i or n_FI must be given a non-empty value');
end

% if both w_b and n_FI are specified, use w_b (throw a warning)
if ~isempty(opts.n_CC) && ~isempty(opts.w_b)
    warning('when both n_CC and w_b are given, n_FI will be ignored');
    opts.n_FI = [];
end

if ~isempty(opts.n_FI)
    % n_CC should be a positive integer
    if ~isa(opts.n_FI,'numeric')
        error('n_FI must be numeric')
    elseif rem(opts.n_FI,1) || opts.n_FI<0
        error('n_FI must be a positive integer');
    end
end

% w_i should be a positive numeric value
if ~isempty(opts.w_i)
    if ~isa(opts.w_i,'numeric')
        error('w_i must be numeric')
    elseif opts.w_i<0
        error('w_i must be positive');
    end
end

%% perform checks for Boundary reduction parameters
% ======================================================================= %

% if Boundary method is "none", w_b will be ignored
if strcmpi(opts.BoundaryMethod,'none') && ~isempty(opts.w_b)
    warning('Boundary reduction type is "none" so w_b will be ignored');
    opts.w_b = [];    
end   
    
% if Boundary method is "none", n_CC will be ignored
if strcmpi(opts.BoundaryMethod,'none') && ~isempty(opts.n_CC)
    warning('Boundary reduction type is "none" so n_CC will be ignored');
    opts.n_CC = [];
end

% neither w_b nor n_CC is given and Boundary method is not "none"
if isempty(opts.w_b) && isempty(opts.n_CC) && ~strcmpi(opts.BoundaryMethod,'none')
    if ~isempty(opts.w_i)
        warning('Either w_b or n_CC should be given a non-empty value unless BoundaryMethod is "none". Setting w_b = w_i and proceeding anyway');
        opts.w_b = opts.w_i;
    else
        error('Either w_b or n_CC must be given a non-empty value unless BoundaryMethod is "none"');
    end
end

if ~isempty(opts.n_CC)
    % n_CC should be a positive integer
    if ~isa(opts.n_CC,'numeric')
        error('n_CC must be numeric')
    elseif rem(opts.n_CC,1) || opts.n_CC<0
        error('n_CC must be a positive integer');
    end
end

% w_b should be a positive numeric value
if ~isempty(opts.w_b)
    if ~isa(opts.w_b,'numeric')
        error('w_b must be numeric')
    elseif opts.w_b<0
        error('w_b must be positive');
    end
end

% if fixed-interface modes are also specified, give a warning
if ~isempty(opts.n_CC) && ~isempty(opts.w_b)
    warning('when both n_CC and w_b are given, n_CC will be ignored');
    opts.n_CC = [];    
end
    
%% additional option checking

% output details of reduction?
if ~isa(opts.verbose,'logical')
    error('verbose must be a logical');
end
    
% plot details of reduction?
if ~isa(opts.plots,'logical')
    error('plots must be a logical');
end