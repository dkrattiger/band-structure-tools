function [K_LIR,M_LIR,dof_sets_LIR,T_LIR] = LIR(K_CB,M_CB,dof_sets,options)

% Dimitri Krattiger
%
% Description
% ===========
% This code uses local interface reduction (LIR) to reduce the interface
% partition of Craig-Bampton mass and stiffness matrices
%
% inputs
% ======
% K_CB      = free Craig-Bampton stiffness matrix
% 
% M_CB      = free Craig-Bampton mass matrix
%
% dof_sets  = DOF set structure for full model
% 
% options   = structure containing options parameters. Some options must
%             be specified:
%
%                1) options.n_CC -or- options.w_b  
%                    (unless options.BoundaryMethod = 'none')
%
%             NEED TO ELABORATE ON OPTIONS 
%
% outputs
% =======
% K_LIR         = free interface-reduced CB stiffness matrix
% 
% M_LIR         = free interface-reduced CB mass matrix
%
% dof_sets_LIR  = updated DOF set structure to go with LIR reduced model
% 
% T_LIR         = transformation between CB DOF vector and CB-LIR DOF vector

%% Check what inputs are given and set rest to default
% ======================================================================= %

% set_default values for options
defaults.BoundaryMethod         = 'exact';
defaults.n_CC                   = [];
defaults.w_b                    = [];
defaults.wQS                    = 0;
defaults.normTypeLIR            = 'mass';
defaults.verbose                = false;
defaults.verboseTab             = '%%   ';
defaults.verboseTabSum          = '';
defaults.plots                  = false;
defaults.outputT                = true;
defaults.preOrthoTypeLIRWeak    = 'none';
% defaults.svdTolLIRExact         = 1e-10;     % This causes a loss of monotonicity but can speed things up

defaults.svdTolLIRWeak          = 1e-3;
defaults.preOrtho               = true;

defaults.orthoTypeLIRExact      = 'svd';     % 'none', 'qr', 'svd'
defaults.svdTolLIRExact         = 0;         % 0 = no reduction

defaults.orthoTypeLIRHybrid     = 'svd';     % 'none', 'qr', 'svd'
defaults.svdTolLIRHybrid        = 1e-6;      % 0 = no reduction

% if doing weak interface method with qr pre-orthogonalization
if exist('options')
    if isfield(options,'BoundaryMethod')
        if isequal(options.BoundaryMethod,'weak')
            if isfield(options,'preOrthoTypeLIRweak')
                if ~isequal(options.preOrthoTypeLIRWeak,'none')                    
                    % no huge difference it seems between these
%                     defaults.svdTolLIRWeak = 1e-3;
%                     defaults.svdTolLIRWeak = 5e-4;
%                     defaults.svdTolLIRWeak = 1e-4;
                    defaults.svdTolLIRWeak = 1e-1;
                end
            end
        end
    end
end


% overwrite defaults with options that are specified
options = setstructfields(defaults,options);

%% Display Timing info
% ======================================================================= %

if options.verbose
    t_start_LIR = tic;
    tb = options.verboseTabSum;
    options.verboseTabSum = [options.verboseTabSum,options.verboseTab];
    tbi = options.verboseTabSum;
    fprintf([tb,'\n'])
    fprintf([tb,'Starting Local Interface Reduction (LIR)\n']) 
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
end

%% check if K_CB is a structure 
% ======================================================================= %
% (if so it indicates that the reduced model has residual 
% enhanced matrices included)
K_CB_full = K_CB;
M_CB_full = M_CB;

res_enhanced = isstruct(K_CB_full);
if res_enhanced
    K_CB = K_CB_full.w0;
    M_CB = M_CB_full.w0;
end

%% extract interior and boundary dof index sets from dof_sets structure
% ======================================================================= %

% take index of all interior DOFs from "dof_sets" structure
i_FI = dof_sets.i';
n_FI = length(i_FI);

% create a list of all boundarey DOFs
dof_set_names = fieldnames(dof_sets);
i_b = [];
for i = 2:length(dof_set_names)
    i_b = [i_b,dof_sets.(dof_set_names{i})'];
end
n_b = length(i_b);

%% Mode Calculation Parameters 
% ======================================================================= %

% if w_cut is available, we ignore n_CC
full_eig = true;
if isempty(options.w_b) | ~isempty(options.n_CC)
    full_eig = false;
elseif options.w_b<0
    full_eig = false;
end

% if 
if ~full_eig & length(i_b)<1500
    full_eig = true;
end

wCenter = options.wQS(ceil(length(options.wQS)/2));

% boundary sets to be paired
boundary_set_names = {{'l','r'},{'f','b'},{'d','t'},...
        {'lf','rf','rb','lb'},{'ld','rd','rt','lt'},{'fd','bd','bt','ft'},...
       {'lfd','rfd','rbd','lbd','lft','rft','rbt','lbt'}};


%% CC Mode Calculation
% ======================================================================= %

if strcmpi(options.BoundaryMethod,'exact') | strcmpi(options.BoundaryMethod,'weak')
    if full_eig
        [PHI_b,L_b] = eig(full(K_CB(i_b,i_b)),...
            full(M_CB(i_b,i_b)));
        [L_b,i_sort] = sort(diag(L_b));    
        PHI_b = PHI_b(:,i_sort);
        if ~isempty(options.n_CC)
            n_CC = options.n_CC;
        else
            n_CC = sum(abs(sqrt(L_b))<=options.w_b);    
        end
    else
        n_CC = options.n_CC;
        [PHI_b,L_b] = eigs(K_CB(i_b,i_b),...
            M_CB(i_b,i_b),n_CC,wCenter^2);
        [L_b,i_sort] = sort(diag(L_b));
        PHI_b = PHI_b(:,i_sort);
    end

    % test for rigid body modes and add them if necessary
    if (min(L_b)>max(L_b)*1e-10) && wCenter == 0;
        PHIr = zeros(length(i_b),3);
        PHIr(1:3:end,1) = 1;
        PHIr(2:3:end,2) = 1;
        PHIr(3:3:end,3) = 1;

        PHI_b = [PHIr,PHI_b];
        L_b = [0;0;0;L_b];
    end


    % sort eigenvalues and eigenvectors by distance from centering frequency
    [~,i_sort] = sort(abs(L_b-wCenter^2));
    L_b = L_b(i_sort);
    PHI_b = PHI_b(:,i_sort);

    % truncate fixed-interface modes and frequencies to n_CC
    L_b = L_b(1:n_CC);
    PHI_b = PHI_b(:,1:n_CC);

    % do a normal sort of fixed-interface modes (ignoring centering freqency)
    [L_b,i_sort] = sort(L_b);
    PHI_b = PHI_b(:,i_sort);

    % Normalize boundary mode shapes
    switch options.normTypeLIR
        case 'mass'
            norm_mat = M_CB(i_b,i_b);
        case 'stiffness'
            norm_mat = K_CB(i_b,i_b);
        case 'unit'
            norm_mat = 1;
    end    
    PHI_b = PHI_b*diag(diag(PHI_b'*norm_mat*PHI_b).^(-1/2));

    if options.preOrtho
        [PHI_b,~] = qr(PHI_b,0);
    end
end

%% Define modes to be used at each interface set
% ======================================================================= %
switch options.BoundaryMethod
    case 'exact'
        % collect boundary mode sets in a structure
        for i = 1:length(boundary_set_names)

            phi_LIR.(boundary_set_names{i}{1}) = [];
            for j = 1:length(boundary_set_names{i})
                phi_LIR.(boundary_set_names{i}{1}) = [phi_LIR.(boundary_set_names{i}{1}),...
                    PHI_b(dof_sets.(boundary_set_names{i}{j})-n_FI,:)];
            end
        end  

        % ensure that modes are well spaced
        switch options.orthoTypeLIRExact
            case 'svd'

                % for each set of boundary modes perform SVD and keep only
                % modes with high enough singular values
                for i = 1:length(boundary_set_names)                
                    [U,S,~] = svd(phi_LIR.(boundary_set_names{i}{1}),'econ');
                    
                    S = diag(S);
                    
                    % define tolerance for cutoff
                    [m,n] = size(phi_LIR.(boundary_set_names{i}{1}));
                    singular_tol = max(m,n) * max(S) * options.svdTolLIRExact; 

                    phi_LIR.(boundary_set_names{i}{1}) = U(:,S > singular_tol);
                end
                
            case 'qr'
                % for each boundary set orthornomormalize the modes using
                % graham-schmidt
                for i = 1:length(boundary_set_names)
                    [phi_LIR.(boundary_set_names{i}{1}),~] = qr(phi_LIR.(boundary_set_names{i}{1}),0); 
                end
            
            case 'none'
                % if too many modes exist for a set to be linearly
                % independent, switch it to identity
                for i = 1:length(boundary_set_names)
                    [m,n] =  size(phi_LIR.(boundary_set_names{i}{1}));
                    if n>m
                    	phi_LIR.(boundary_set_names{i}{1}) = eye(m);
                    end
                end
                warning('Ill conditioning may cause bad solutions when no LIR modes are not orthonormalized')
        end

        % find size of boundary mode sets
        n_master = zeros(1,length(boundary_set_names));
        for i = 1:length(boundary_set_names)
            n_master(i) = size(phi_LIR.(boundary_set_names{i}{1}),2);
        end
        
    case 'weak'
        
        % collect boundary mode sets in a structure
        n_master = zeros(1,length(boundary_set_names));
        
        for i = 1:length(boundary_set_names)
            % number of DOFs in current boundary set
            n_dof_set = length(dof_sets.(boundary_set_names{i}{1}));
            
            % number of connecting boundary sets
            n_sets = length(boundary_set_names{i});           
            
            if n_dof_set<n_CC
                n_master(i) = n_dof_set;
                for j = 1:n_sets
                    phi_LIR.(boundary_set_names{i}{j}) = eye(n_dof_set);
                end
            else
                
                % total number of modal dofs
                n_mode = n_sets*n_CC;

                % total number of physical constraints
                n_con = (n_sets-1)*n_dof_set;

                % preallocate constraint matrix
                a = zeros(n_con,n_mode);
                row_count = 0;
                col_count = 0;

                % orthogonalize boundary modes
                PHI_master = PHI_b(dof_sets.(boundary_set_names{i}{1})-n_FI,1:n_CC);
                switch options.preOrthoTypeLIRWeak
                    case 'svd'
                        [PHI_master,~,~] = svd(PHI_master,0);
                    case 'qr'                        
                        [PHI_master,~] = qr(PHI_master,0);
                    case 'none'
                end
                PHI_b(dof_sets.(boundary_set_names{i}{1})-n_FI,1:n_CC) = PHI_master;
                
                % place master modes into constraint matrix
                a(1:n_con,(1:n_CC)) = repmat(PHI_master,[(n_sets-1),1]);

                col_count = col_count + n_CC;

                for j = 2:n_sets
                    PHI_slave = PHI_b(dof_sets.(boundary_set_names{i}{j})-n_FI,1:n_CC);
                    switch options.preOrthoTypeLIRWeak
                        case 'svd'
                            [PHI_slave,~,~] = svd(PHI_slave,0);
                        case 'qr'
                            [PHI_slave,~] = qr(PHI_slave,0);
                        case 'none'
                    end
                    PHI_b(dof_sets.(boundary_set_names{i}{j})-n_FI,1:n_CC) = PHI_slave;
                    a(row_count+(1:n_dof_set),col_count+(1:n_CC)) = -PHI_slave;

                    col_count = col_count + n_CC;
                    row_count = row_count + n_dof_set;
                end

                % compute svd of constraint matrix                
                [~,Sa,Va] = svd(a,0);
                Sa = [diag(Sa)];

                % define tolerance and use it to select "master" generalized
                % coordinates
                [m,n] = size(a);
                %singular_tol = max(m,n) * max(Sa) * options.svdTolLIRWeak; 
                singular_tol = max(Sa) * options.svdTolLIRWeak; 

                % number of master gen. coords
                n_slave(i) = sum(Sa>=singular_tol);
                n_master(i) = n_mode-n_slave(i);         
                %n_master(i) = min(sum(Sa<tol),n_CC);   

                if n_master(i) >= n_dof_set
                    for j = 1:n_sets
                        phi_LIR.(boundary_set_names{i}{j}) = eye(n_dof_set);
                    end
                    n_master(i) = n_dof_set;
                else

                    %i_master = (n_mode-n_LI+1):n_mode;
                    i_master = (n_slave(i)+1):n_mode;
                    
                    B = Va(:,i_master);
                    for j = 1:n_sets
                        phi_LIR.(boundary_set_names{i}{j}) = ...
                            PHI_b(dof_sets.(boundary_set_names{i}{j})-n_FI,:)*...
                            B(((j-1)*n_CC+1):(j*n_CC),:);
                    end
                end
            end
        end
        
    case 'hybrid'
        
        % create DOF sets structure that is empty everywhere
        dof_sets_clear = dof_sets;
        dof_set_names = fieldnames(dof_sets_clear);
        for j = 1:length(dof_set_names)
            dof_sets_clear.(dof_set_names{j}) = [];
        end
        
        % partially couple model boundary mode sets
        for i = 1:length(boundary_set_names)
            
            if options.n_CC >= length(dof_sets.(boundary_set_names{i}{1}))
                phi_LIR.(boundary_set_names{i}{1}) = ...
                    eye(length(dof_sets.(boundary_set_names{i}{1})));
            else
                dof_sets_temp = dof_sets_clear;
                %dof_sets_temp.i = (1:n_FI)';
                
                % If DOFs are in active boundary set, add them into temporary
                % DOF set structure. Otherwise add them into interior field.
                for j = 1:length(boundary_set_names)
                    if i==j
                        for k = 1:length(boundary_set_names{j})
                            dof_sets_temp.(boundary_set_names{j}{k}) = ...
                                 dof_sets.(boundary_set_names{j}{k});
                        end
                    else
                        for k = 1:length(boundary_set_names{j})
                            dof_sets_temp.i = [dof_sets_temp.i;...
                                 dof_sets.(boundary_set_names{j}{k})];
                        end
                    end
                end

                % note that the FI modes are not added into the dof_sets_temp
                % structure

                % form transformation to couple just the active boundary set
                T_per = Periodic_Boundary_Conditions(dof_sets_temp);
                T_per_gamma = T_per.s0 + ...
                              T_per.s1 + T_per.s2 + T_per.s3 + ...
                              T_per.s12 + T_per.s13 + T_per.s23 + T_per.s123;

                % apply periodicity transform (gamma because zero phase offset)
                Kb = T_per_gamma'*K_CB((n_FI+1):end,(n_FI+1):end)*T_per_gamma;
                Mb = T_per_gamma'*M_CB((n_FI+1):end,(n_FI+1):end)*T_per_gamma;
                
                % ensure symmetry
                Kb = (1/2)*(Kb+Kb');
                Mb = (1/2)*(Mb+Mb');
                
%                 doHybridCondensation = false;
%                 if doHybridCondensation
%                     % perform static condensation of unnecessary boundary DOFs
%                     % onto master boundary DOFs
%                     i_s = 1:length(dof_sets_temp.i);
%                     i_m = (length(dof_sets_temp.i)+1):size(Kb,1);
%                     Tcond = [-Kb(i_s,i_s)\Kb(i_s,i_m);...
%                              eye(length(i_m))];
% 
%                     % % apply static condensation to boundary segment mass and
%                     % stiffness matrices
%                     Kb = Tcond'*Kb*Tcond;
%                     Mb = Tcond'*Mb*Tcond;
%                 end
                
                % compute eigenvalues
                if full_eig
                    [PHI_b,L_b] = eig(full(Kb),full(Mb));
                    L_b = diag(L_b);
                    
                    if isempty(options.n_CC)
                        n_CC = sum(abs(sqrt(L_b))<=options.w_b);
                    else
                        n_CC = options.n_CC;
                    end
                    
                else
                    n_CC = options.n_CC;
                    [PHI_b,L_b] = eigs(Kb,Mb,n_CC,wCenter^2);
                    L_b = diag(L_b);
                end


                % % sort eigenvalues and eigenvectors by distance from centering frequency
                [~,i_sort] = sort(abs(L_b-wCenter^2));
                L_b = L_b(i_sort);
                PHI_b = PHI_b(:,i_sort);
                
                % truncate fixed-interface modes and frequencies to n_CC
                L_b = L_b(1:n_CC);
                PHI_b = PHI_b(:,1:n_CC);

                % do a normal sort of fixed-interface modes (ignoring centering freqency)
                [L_b,i_sort] = sort(L_b);
                PHI_b = PHI_b(:,i_sort);
%                 
%                 L_b
%                 pause

                % Normalize boundary mode shapes
                switch options.normTypeLIR
                    case 'mass'
                        norm_mat = Mb;
                    case 'stiffness'
                        norm_mat = Kb;
                    case 'unit'
                        norm_mat = 1;
                end
                PHI_b = PHI_b*diag(diag(PHI_b'*norm_mat*PHI_b).^(-1/2));


                % by the way Kb and Mb are defined, the interesting (useful) part of
                % the mode shape is always contained right at the end of PHI_b
%                 if doHybridCondensation
%                     phi_LIR.(boundary_set_names{i}{1}) = PHI_b;
%                 else
                phi_LIR.(boundary_set_names{i}{1}) = PHI_b((length(dof_sets_temp.i)+1):size(PHI_b,1),:);
%                 end
                
                if false
                    boundary_set_names{i}{1}
                    figure(99);clf
                    PHIplot = phi_LIR.(boundary_set_names{i}{1});
                    xcoords = 1:size(PHIplot,1)/2;
                    ycoords = 1:size(PHIplot,2);
                    [xgrid,ygrid] = ndgrid(xcoords,ycoords);

                    PHIplot = PHIplot*diag(diag(PHIplot'*PHIplot).^(-1/2))*1.5;

                    plot(ygrid,xgrid,'k.-');hold on
                    plot(ygrid+PHIplot(2:2:end,:),xgrid+PHIplot(1:2:end,:),'o--')
                    pause
                end
                
                % ensure that modes are well spaced
                switch options.orthoTypeLIRHybrid
                    case 'svd'
                        
                        % for each set of boundary modes perform SVD and keep only
                        % modes with high enough singular values  
                        [U,S,~] = svd(phi_LIR.(boundary_set_names{i}{1}),'econ');

                        S = diag(S);
                        

                        % define tolerance for cutoff
                        [m,n] = size(phi_LIR.(boundary_set_names{i}{1}));
                        singular_tol = max(m,n) * max(S) * options.svdTolLIRHybrid; 

                        phi_LIR.(boundary_set_names{i}{1}) = U(:,S > singular_tol);

%                         figure(888);clf;
%                         semilogy(S);hold on
%                         semilogy(S(S > singular_tol),'r.')
%                         sum((S > singular_tol))
%                         pause
                        
                    case 'qr'
                        % for each boundary set orthornomormalize the modes using
                        % graham-schmidt
                        [phi_LIR.(boundary_set_names{i}{1}),~] = qr(phi_LIR.(boundary_set_names{i}{1}),0); 

                    case 'none'
                        % if too many modes exist for a set to be linearly
                        % independent, switch it to identity
                        [m,n] =  size(phi_LIR.(boundary_set_names{i}{1}));
                        if n>m
                            phi_LIR.(boundary_set_names{i}{1}) = eye(m);
                        end
                        warning('Ill conditioning may cause bad solutions when no LIR modes are not orthonormalized')
                end
                
            end
        end
        
        % find size of boundary mode sets
        n_master = zeros(1,length(boundary_set_names));
        for i = 1:length(boundary_set_names)
            n_master(i) = size(phi_LIR.(boundary_set_names{i}{1}),2);
        end
        
end

%% Form Local interface reduction tranformation
% ======================================================================= %

% sum up interface DOFs
n_b_LIR = n_master(1)*2 + n_master(2)*2 + n_master(3)*2 + ...
          n_master(4)*4 + n_master(5)*4 + n_master(6)*4 + ...
          n_master(7)*8;
n_b_LIR_per = sum(n_master);
n_dof_LIR = n_FI + n_b_LIR;

n_dof_BMS = n_b + n_FI;
T_LIR = zeros(n_dof_BMS,n_dof_LIR);

% T_LIR(1:n_FI,1:n_FI) = eye(n_FI);
dof_sets_LIR = dof_sets;
count = n_FI;
for i = 1:length(boundary_set_names)
    for j = 1:length(boundary_set_names{i})
        dof_sets_LIR.(boundary_set_names{i}{j}) = (count+1:count+n_master(i))';
        count = count+n_master(i);
    end
end  

T_LIR(dof_sets.i,dof_sets_LIR.i) = eye(n_FI);
for i = 1:length(boundary_set_names)
    for j = 1:length(boundary_set_names{i})
        switch options.BoundaryMethod
            case 'exact'
                ind = 1;
            case 'weak'
                ind = j;
            case 'hybrid'
                ind = 1;
        end
        T_LIR(dof_sets.(boundary_set_names{i}{j}),...
            dof_sets_LIR.(boundary_set_names{i}{j})) = ...
            phi_LIR.(boundary_set_names{i}{ind});
    end
end  

T_LIR = sparse(T_LIR);
if res_enhanced
    M_LIR.w0 = T_LIR'*M_CB_full.w0*T_LIR;
    K_LIR.w0 = T_LIR'*K_CB_full.w0*T_LIR;
    M_LIR.w2 = T_LIR'*M_CB_full.w2*T_LIR;
    K_LIR.w2 = T_LIR'*K_CB_full.w2*T_LIR;
    M_LIR.w4 = T_LIR'*M_CB_full.w4*T_LIR;
    K_LIR.w4 = T_LIR'*K_CB_full.w4*T_LIR;
else
    M_LIR = T_LIR'*M_CB*T_LIR;
    K_LIR = T_LIR'*K_CB*T_LIR;
end

if ~options.outputT
    T_LIR = [];
end

%% Display Timing info
% ======================================================================= %

if options.verbose
    fprintf([tbi,'LIR calculation time: %5.2f s\n'],toc(t_start_LIR))
    fprintf([tb,repmat('%%',1,60-length(sprintf(tb))),'\n'])
    fprintf([tb,'\n'])
end