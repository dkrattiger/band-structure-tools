function [node_sets] = find_node_sets(X,R)
%
% Dimitri Krattiger (4-3-2017)
%
% description
% ===========
% This code computes node sets structure by identifying nodes that are
% separated exactly by one lattice vector
%
% inputs
% ======
%
% X         = nodal coordinates for FEM model 
%            (X=[x(:)] in 1D models
%             X=[x(:),y(:)] in 2D models
%             X=[x(:),y(:),z(:)] in 3D models)
%
% R         = lattice vectors in 1, 2, or 3 directions 
%             (R = r1 for 1 direction of periodicity i.e. 1P
%              R = [r1,r2] for 2 directions of periodicity i.e. 2P
%              R = [r1,r2,r3] for 3 directions of periodicity i.e. 3P)
%
% outputs
% =======
% node_sets = structure containing node sets


%% Represent Coordinates in lattice-vector coordinates
% ======================================================================= %

% number of significant digits
% sig_digits = 4;
% round_val = 10^(sig_digits-ceil(log10(max(abs(X(:))))));


% get size of coordinate array and transpose if necessary
[a,b] = size(X);
if a<b
    X = X';
end
n_p = size(R,2);
[n_nodes,n_dim_x] = size(X);

% Normalize X and R
R = R/max(abs(X(:)));
X = X/max(abs(X(:)));

% round X to make identifying easier
% X = round(X*round_val)/round_val;
% R = round(R*round_val)/round_val;

% [Lia,Locb] = ismember(X+ones(n_nodes,1)*R(:,1)',X,'rows');
tolmult = 1e-6;
[Lia,Locb] = ismember_tol(X+ones(n_nodes,1)*R(:,1)',X,'rows','tolerance',tolmult);
i_l1 = find(Lia);
i_r1 = Locb(Lia);


% (Lia)
% (Locb)

% i_l1
% i_r1
% pause

% figure(99)
% plot3(X(:,1),X(:,2),X(:,3),'k.');hold on
% plot3(X(i_l1,1),X(i_l1,2),X(i_l1,3),'ro');
% plot3(X(i_r1,1),X(i_r1,2),X(i_r1,3),'go');
% 
% size(i_l1)
% size(i_r1)
% pause

% [~,i_sort] = ismember(X(i_l1,:)+ones(length(i_l1),1)*R(:,1)',X(i_r1,:),'rows');
[~,i_sort] = ismember_tol(X(i_l1,:)+ones(length(i_l1),1)*R(:,1)',X(i_r1,:),'rows','tolerance',tolmult);

%i_r1 = i_r1(i_sort);

i_f1 = [];  i_b1 = [];  i_d1 = [];  i_t1 = [];
if n_p >= 2
%     [Lia,Locb] = ismember(X+ones(n_nodes,1)*R(:,2)',X,'rows');
    [Lia,Locb] = ismember_tol(X+ones(n_nodes,1)*R(:,2)',X,'rows','tolerance',tolmult);
    i_f1 = find(Lia);
    i_b1 = Locb(Lia);
%     [~,i_sort] = ismember(X(i_f1,:)+ones(length(i_f1),1)*R(:,2)',X(i_b1,:),'rows');
    [~,i_sort] = ismember_tol(X(i_f1,:)+ones(length(i_f1),1)*R(:,2)',X(i_b1,:),'rows','tolerance',tolmult);
    i_b1 = i_b1(i_sort);
end

if n_p >= 3
%     [Lia,Locb] = ismember(X+ones(n_nodes,1)*R(:,3)',X,'rows');
    [Lia,Locb] = ismember_tol(X+ones(n_nodes,1)*R(:,3)',X,'rows','tolerance',tolmult);
    i_d1 = find(Lia);
    i_t1 = Locb(Lia);
%     [~,i_sort] = ismember(X(i_d1,:)+ones(length(i_d1),1)*R(:,3)',X(i_t1,:),'rows');
    [~,i_sort] = ismember_tol(X(i_d1,:)+ones(length(i_d1),1)*R(:,3)',X(i_t1,:),'rows','tolerance',tolmult);
    i_t1 = i_t1(i_sort);
end

% arrange each
i_l1 = i_l1(:)';
i_r1 = i_r1(:)';
i_f1 = i_f1(:)';
i_b1 = i_b1(:)';
i_d1 = i_d1(:)';
i_t1 = i_t1(:)';

% [i_l1;i_r1]
% [i_f1;i_b1] 
% [i_d1;i_t1] 
% pause

%% determine boundary sets by analyzing nodal coordinates
% ======================================================================= %

% interior
i_i = setdiff(1:n_nodes,[i_l1,i_r1,i_f1,i_b1,i_d1,i_t1],'stable');

% left and right faces
i_l = setdiff(i_l1,[i_f1,i_b1,i_d1,i_t1],'stable');
i_r = setdiff(i_r1,[i_f1,i_b1,i_d1,i_t1],'stable');

% front and back faces
i_f = setdiff(i_f1, [i_l1,i_r1,i_d1,i_t1],'stable');
i_b = setdiff(i_b1, [i_l1,i_r1,i_d1,i_t1],'stable');

% down and top faces
i_d = setdiff(i_d1, [i_l1,i_r1,i_f1,i_b1],'stable');
i_t = setdiff(i_t1, [i_l1,i_r1,i_f1,i_b1],'stable');

% left-front, right-front, right-back and left-back edges
i_lf = setdiff(intersect(i_l1,i_f1,'stable'),[i_d1,i_t1],'stable');
i_rf = setdiff(intersect(i_r1,i_f1,'stable'),[i_d1,i_t1],'stable');
i_rb = setdiff(intersect(i_r1,i_b1,'stable'),[i_d1,i_t1],'stable');
i_lb = setdiff(intersect(i_l1,i_b1,'stable'),[i_d1,i_t1],'stable');

% left-down, right-down, right-top, and left-top edges
i_ld = setdiff(intersect(i_l1,i_d1,'stable'),[i_f1,i_b1],'stable');
i_rd = setdiff(intersect(i_r1,i_d1,'stable'),[i_f1,i_b1],'stable');
i_rt = setdiff(intersect(i_r1,i_t1,'stable'),[i_f1,i_b1],'stable');
i_lt = setdiff(intersect(i_l1,i_t1,'stable'),[i_f1,i_b1],'stable');

% front-down, back-down, back-top, and front-top edges
i_fd = setdiff(intersect(i_f1,i_d1,'stable'),[i_l1,i_r1],'stable');
i_bd = setdiff(intersect(i_b1,i_d1,'stable'),[i_l1,i_r1],'stable');
i_bt = setdiff(intersect(i_b1,i_t1,'stable'),[i_l1,i_r1],'stable');
i_ft = setdiff(intersect(i_f1,i_t1,'stable'),[i_l1,i_r1],'stable');

% left-front-down, right-front-down, right-back-down, and left-back-down
% corners
i_lfd = intersect(intersect(i_l1,i_f1,'stable'),i_d1,'stable');
i_rfd = intersect(intersect(i_r1,i_f1,'stable'),i_d1,'stable');
i_rbd = intersect(intersect(i_r1,i_b1,'stable'),i_d1,'stable');
i_lbd = intersect(intersect(i_l1,i_b1,'stable'),i_d1,'stable');

% left-front-top, right-front-top, right-back-top, and left-back-top
% corners
i_lft = intersect(intersect(i_l1,i_f1,'stable'),i_t1,'stable');
i_rft = intersect(intersect(i_r1,i_f1,'stable'),i_t1,'stable');
i_rbt = intersect(intersect(i_r1,i_b1,'stable'),i_t1,'stable');
i_lbt = intersect(intersect(i_l1,i_b1,'stable'),i_t1,'stable');

%% collect node index sets into structure
% ======================================================================= %

% interior
node_sets.i = i_i;

% left and right faces
node_sets.l = i_l;
node_sets.r = i_r;

% front and back faces
node_sets.f = i_f;
node_sets.b = i_b;

% down and top faces
node_sets.d = i_d;
node_sets.t = i_t;

% left-front, right-front, right-back and left-back edges
node_sets.lf = i_lf;
node_sets.rf = i_rf;
node_sets.rb = i_rb;
node_sets.lb = i_lb;

% left-down, right-down, right-top, and left-top edges
node_sets.ld = i_ld;
node_sets.rd = i_rd;
node_sets.rt = i_rt;
node_sets.lt = i_lt;

% front-down, back-down, back-top, and front-top edges
node_sets.fd = i_fd;
node_sets.bd = i_bd;
node_sets.bt = i_bt;
node_sets.ft = i_ft;

% left-front-down, right-front-down, right-back-down, and left-back-down
% corners
node_sets.lfd = i_lfd;
node_sets.rfd = i_rfd;
node_sets.rbd = i_rbd;
node_sets.lbd = i_lbd;

% left-front-top, right-front-top, right-back-top, and left-back-top
% corners
node_sets.lft = i_lft;
node_sets.rft = i_rft;
node_sets.rbt = i_rbt;
node_sets.lbt = i_lbt;


function [Lia,Lib] = ismember_tol(A,B,varargin)

% ismember function that determines set membership but within some
% tolerance

option_args = varargin;
for j = (nargin-3):-1:1
    if strcmpi(varargin{j},'tolerance')
        tol = option_args{j+1};
        option_args([j,j+1]) = [];
    end
end

% rounding digits;
rd = 10^(-ceil(log10(tol)));

[Lia,Lib] = ismember(round(A*rd)/rd,round(B*rd)/rd,option_args{:});