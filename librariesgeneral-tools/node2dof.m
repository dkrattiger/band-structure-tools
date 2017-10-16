function output = node2dof(input,varargin)
%
% Dimitri Krattiger (4-3-2017)
%
% description
% ===========
% This code converts node indices to DOF indices and supports various
% input/output types.
%
% inputs
% ======
%
% input     = either (1) a list of node indices, or 
%                    (2) a structure whose entries are lists of node
%                        indices
%
% n_dpn     = scalar: number of DOFs per node for all nodes
%              -or-
%             vector: where ith entry gives number of DOFs for ith node
%              -or-
% DOF       = vector which for every node contains (node#).(DOF#)
%             (e.g. if node 1234 has 6 DOFs:
%              DOF = [...
%                     1234.1;   (x-displacement of node 1234)
%                     1234.2;   (y-displacement of node 1234)
%                     1234.3;   (z-displacement of node 1234)
%                     1234.4;   (x-rotation of node 1234)
%                     1234.5;   (y-rotation of node 1234)
%                     1234.6;   (z-rotation of node 1234)
%                     ...]
%
%
% outputs
% =======
% output    = either (1) a list of dof indices, or 
%                    (2) a structure whose entries are lists of dof
%                        indices       
        
if isstruct(input)
    
    % if input is a structure, use recursive calls to "node2dof" to treat
    % each individual field of the structure
    node_sets = input;
    node_set_names = fieldnames(node_sets);
        
    % loop through field names and recursively apply node2dof
    for i = 1:length(node_set_names)
        dof_sets.(node_set_names{i}) = ...
            node2dof(node_sets.(node_set_names{i}),varargin{:});
    end
        
    output = dof_sets;
    
else
    
    i_node = input(:);
    n_node = length(i_node);

    if ~isempty(i_node)
        if nargin == 2
            n_dpn = varargin{1};

            if length(n_dpn) == 1

                % allocate space for dof indices
                i_dof = zeros(n_dpn,length(i_node));

                % quickest way to assign many dof indices at once
                for i = 1:n_dpn
                    if ~isempty(i_node)
                        i_dof(i,:) = (n_dpn)*(i_node-1) + i;
                    end
                end
                i_dof = reshape(i_dof,[numel(i_dof),1]);
            else  

                % cumulative sum of DoFs-per-node vector
                cumsum_n_dpn = [0,cumsum(n_dpn)];

                % determine size of new DOF index
                n_dof = sum(n_dpn(i_node));

                % initialize and preallocate
                i_dof = zeros(n_dof,1);
                count_DOF = 0;
                for i = 1:n_node

                    % sum up number of dofs for nodes below current node
                    count_prev = cumsum_n_dpn(i_node(i));

                    % assign DOF indices
                    i_dof(count_DOF+(1:n_dpn(i)),:) = count_prev+(1:n_dpn(i));

                    % increment counter
                    count_DOF = count_DOF+n_dpn(i);
                end
            end
        elseif nargin == 3

            % node_names	is a list of all node identifiers in system (all
            %               integers
            %
            % DOF        	is a list of all DOF names in system. The integer
            %               component correponds to that DOF's node name
            %
            % i_node        is the node index that is being converted to a
            %               dof index

            node_names = varargin{1};
            DOF = varargin{2};

            % guess for maximum size of i_dof in order to preallocate
            i_dof = zeros(length(i_node)*6*2,1);
            count = 0;
            for i = 1:length(i_node)
                i_dof_add = find(floor(DOF) == node_names(i_node(i)));
                n_dof_add = length(i_dof_add);
                i_dof(count+(1:n_dof_add)) = i_dof_add;
                count = count+n_dof_add;            
            end
            
            % truncate extra zeros
            i_dof = i_dof(1:count);
        end
    else
        i_dof = [];
    end
    output = i_dof;
        
end
