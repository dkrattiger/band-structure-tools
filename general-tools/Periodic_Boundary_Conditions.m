function T_per = Periodic_Boundary_Conditions(dof_sets)

%% Documentation
% ======================================================================= %

% Description
% ===========
% This code generates a set of transformation matrices that can be applied
% to a free model to enforce periodic boundary conditions for any given
% wave vector. This code is aimed at 3D models, but will work for 2D models
% if the correct degree of freedom (dof) sets are given.

% subscript notation guide
% ========================
% i = internal
% l = left (x=0)
% r = right (x=Lx)
% f = front (y=0)
% b = back (y=Ly)
% d = bottom (down) (z=0)
% t = top (z=Lz)

% inputs
% ======
% The following lists give the degree of freedom (dof) indices in the
% system. The sets that are paired to be connected (i.e. i_l and i_r)
% should sort the corresponding degrees of freedom so that they match up.
% If the model contains no degrees of freedom for a certain face, edge, or
% corner, then an empty set [] may be given as input
%
% i_i = list of internal dofs
% i_l = list of left face dofs (not including corners or edges)
% i_r = list of right face dofs (not including corners or edges)
% i_f = list of front face dofs (not including corners or edges)
% i_b = list of back face dofs (not including corners or edges)
% i_d = list of bottom face dofs (not including corners or edges)
% i_t = list of top face dofs (not including corners or edges)
% i_lf = list of left-front edge dofs (not including corners)
% i_rf = list of right-front edge dofs (not including corners)
% i_rb = list of right-back edge dofs (not including corners)
% i_lb = list of left-back edge dofs (not including corners)
% i_ld = list of left-bottom edge dofs (not including corners)
% i_rd = list of right-bottom edge dofs (not including corners)
% i_rt = list of right-top edge dofs (not including corners)
% i_lt = list of left-top edge dofs (not including corners)
% i_fd = list of front-bottom edge dofs (not including corners)
% i_bd = list of back-down edge dofs (not including corners)
% i_bt = list of back-top edge dofs (not including corners)
% i_ft = list of front-top edge dofs (not including corners)
% i_lfd = list of left-front-bottom corner dofs
% i_rfd = list of right-front-bottom corner dofs
% i_rbd = list of right-back-bottom corner dofs
% i_lbd = list of left-back-bottom corner dofs
% i_lft = list of left-front-top corner dofs
% i_rft = list of right-front-top corner dofs
% i_rbt = list of right-back-top corner dofs
% i_lbt = list of left-back-top corner dofs

% outputs
% =======
% The outputs are a set of transformation matrices that can be used to
% create an overall transformation for any desired wave-vector.
%
% T_per = T_per.s0 + T_per.s1*lamx + T_per.s2*lamy + T_per.s3*lamz +
%         T_per.s12*lamx*lamy + T_per.s13*lamx*lamz + T_per.s23*lamy*lamz +
%         T_per.s123*lamx*lamy*lamz
% 
% where lamx = exp(i*kx*Lx), lamy = exp(i*ky*Ly), lamz = exp(i*kz*Lz)
%
% The final periodicity transformation T_per can be used to enforce
% periodicity by applying to the free mass and stiffness matrices as
% follows,
%
% M_per = T_per'*M*T_per,   K_per = T_per'*K*T_per,
%
% where ' indicates the complex conjugate transpose. Similarly, the
% periodic dynamical matrix can be formed as follows.
%
% D_per = T_per'*D*T_per
      
%% Begin Code
% ======================================================================= %  
% extract dof sets from dof_set structure
i_i = dof_sets.i(:);
i_l = dof_sets.l(:);
i_r = dof_sets.r(:);
i_f = dof_sets.f(:);
i_b = dof_sets.b(:);
i_d = dof_sets.d(:);
i_t = dof_sets.t(:);

i_lf = dof_sets.lf(:);
i_rf = dof_sets.rf(:);
i_rb = dof_sets.rb(:);
i_lb = dof_sets.lb(:);

i_ld = dof_sets.ld(:);
i_rd = dof_sets.rd(:);
i_rt = dof_sets.rt(:);
i_lt = dof_sets.lt(:);

i_fd = dof_sets.fd(:);
i_bd = dof_sets.bd(:);
i_bt = dof_sets.bt(:);
i_ft = dof_sets.ft(:);

i_lfd = dof_sets.lfd(:);
i_rfd = dof_sets.rfd(:);
i_rbd = dof_sets.rbd(:);
i_lbd = dof_sets.lbd(:);
i_lft = dof_sets.lft(:);
i_rft = dof_sets.rft(:);
i_rbt = dof_sets.rbt(:);
i_lbt = dof_sets.lbt(:);


% number of interior dofs
n_i = length(i_i);

% number of face dofs
n_l = length(i_l);
n_r = length(i_r);
n_f = length(i_f);
n_b = length(i_b);
n_t = length(i_t);
n_d = length(i_d);

% number of edge dofs
n_lf = length(i_lf);
n_lb = length(i_lb);
n_ld = length(i_ld);
n_lt = length(i_lt);
n_rf = length(i_rf);
n_rb = length(i_rb);
n_rd = length(i_rd);
n_rt = length(i_rt);
n_fd = length(i_fd);
n_ft = length(i_ft);
n_bd = length(i_bd);
n_bt = length(i_bt);

% number of corner dofs
n_lfd = length(i_lfd);
n_lft = length(i_lft);
n_lbd = length(i_lbd);
n_lbt = length(i_lbt);
n_rfd = length(i_rfd);
n_rft = length(i_rft);
n_rbd = length(i_rbd);
n_rbt = length(i_rbt);

i_all = [i_i;...
         i_l;i_r;...
         i_f;i_b;...
         i_d;i_t;...
         i_lf;i_rf;i_rb;i_lb;...
         i_ld;i_rd;i_rt;i_lt;...
         i_fd;i_bd;i_bt;i_ft;...
         i_lfd;i_rfd;i_rbd;i_lbd;i_lft;i_rft;i_rbt;i_lbt];

% create index to "un-sort" the transfer matrix rows to match the original
% dof sort
[~,i_unsort1] = sort(i_all);         

% create index to "un-sort" the transfer matrix columns to match the
% original dof sort
i_per = [i_i;i_l;i_f;i_d;i_lf;i_ld;i_fd;i_lfd];
[~,i_unsort2] = sort(i_per);  


n_all = length(i_all);
n_per = length(i_per);

% preallocate arrays
T_per.s0 = sparse(n_all,n_per);
T_per.s1 = sparse(n_all,n_per);
T_per.s2 = sparse(n_all,n_per);
T_per.s3 = sparse(n_all,n_per);
T_per.s12 = sparse(n_all,n_per);
T_per.s13 = sparse(n_all,n_per);
T_per.s23 = sparse(n_all,n_per);
T_per.s123 = sparse(n_all,n_per);

%% map interior to interior
% ======================================================================= %
T_per.s0(1:n_i,1:n_i) = speye(n_i);


%% map faces to corresponding faces
% ======================================================================= %
T_per.s0(n_i+1:n_i+n_l,n_i+1:n_i+n_l) = speye(n_l);
T_per.s1(n_i+n_l+1:n_i+n_l+n_r,n_i+1:n_i+n_l) = speye(n_r);

T_per.s0(n_i+n_l+n_r+1:n_i+n_l+n_r+n_f,n_i+n_l+1:n_i+n_l+n_f) = speye(n_f);
T_per.s2(n_i+n_l+n_r+n_f+1:n_i+n_l+n_r+n_f+n_b,n_i+n_l+1:n_i+n_l+n_f) = speye(n_b);

T_per.s0(n_i+n_l+n_r+n_f+n_b+1:n_i+n_l+n_r+n_f+n_b+n_d,n_i+n_l+n_f+1:n_i+n_l+n_f+n_d) = speye(n_d);
T_per.s3(n_i+n_l+n_r+n_f+n_b+n_d+1:n_i+n_l+n_r+n_f+n_b+n_d+n_t,n_i+n_l+n_f+1:n_i+n_l+n_f+n_d) = speye(n_t);

% intermediate summations
n_int1 = n_i+n_l+n_r+n_f+n_b+n_d+n_t;
n_int2 = n_i+n_l+n_f+n_d;

%% map edges to corresponding edges
% ======================================================================= %
T_per.s0(n_int1+1:n_int1+n_lf,n_int2+1:n_int2+n_lf) = speye(n_lf);
T_per.s1(n_int1+n_lf+1:n_int1+n_lf+n_rf,n_int2+1:n_int2+n_lf) = speye(n_rf);
T_per.s12(n_int1+n_lf+n_rf+1:n_int1+n_lf+n_rf+n_rb,n_int2+1:n_int2+n_lf) = speye(n_rb);
T_per.s2(n_int1+n_lf+n_rf+n_rb+1:n_int1+n_lf+n_rf+n_rb+n_lb,n_int2+1:n_int2+n_lf) = speye(n_lb);

% intermediate summations
n_int3 = n_i+n_l+n_r+n_f+n_b+n_d+n_t+n_lf+n_rf+n_rb+n_lb;
n_int4 = n_i+n_l+n_f+n_d+n_lf;

T_per.s0(n_int3+1:n_int3+n_ld,n_int4+1:n_int4+n_ld) = speye(n_ld);
T_per.s1(n_int3+n_ld+1:n_int3+n_ld+n_rd,n_int4+1:n_int4+n_ld) = speye(n_rd);
T_per.s13(n_int3+n_ld+n_rd+1:n_int3+n_ld+n_rd+n_rt,n_int4+1:n_int4+n_ld) = speye(n_rt);
T_per.s3(n_int3+n_ld+n_rd+n_rt+1:n_int3+n_ld+n_rd+n_rt+n_lt,n_int4+1:n_int4+n_ld) = speye(n_lt);

% intermediate summations
n_int5 = n_i+n_l+n_r+n_f+n_b+n_d+n_t+n_lf+n_rf+n_rb+n_lb+n_ld+n_rd+n_rt+n_lt;
n_int6 = n_i+n_l+n_f+n_d+n_lf+n_ld;

T_per.s0(n_int5+1:n_int5+n_fd,n_int6+1:n_int6+n_fd) = speye(n_fd);
T_per.s2(n_int5+n_fd+1:n_int5+n_fd+n_bd,n_int6+1:n_int6+n_fd) = speye(n_bd);
T_per.s23(n_int5+n_fd+n_bd+1:n_int5+n_fd+n_bd+n_bt,n_int6+1:n_int6+n_fd) = speye(n_bt);
T_per.s3(n_int5+n_fd+n_bd+n_bt+1:n_int5+n_fd+n_bd+n_bt+n_ft,n_int6+1:n_int6+n_fd) = speye(n_ft);

% intermediate summations
n_int7 = n_i+n_l+n_r+n_f+n_b+n_d+n_t+n_lf+n_rf+n_rb+n_lb+n_ld+n_rd+n_rt+n_lt+n_fd+n_bd+n_bt+n_ft;
n_int8 = n_i+n_l+n_f+n_d+n_lf+n_ld+n_fd;

%% map corners to corresponding corners
% ======================================================================= %
T_per.s0(n_int7+1:n_int7+n_lfd,n_int8+1:n_int8+n_lfd) = speye(n_lfd);
T_per.s1(n_int7+n_lfd+1:n_int7+n_lfd+n_rfd,n_int8+1:n_int8+n_lfd) = speye(n_rfd);
T_per.s12(n_int7+n_lfd+n_rfd+1:n_int7+n_lfd+n_rfd+n_rbd,n_int8+1:n_int8+n_lfd) = speye(n_rbd);
T_per.s2(n_int7+n_lfd+n_rfd+n_rbd+1:n_int7+n_lfd+n_rfd+n_rbd+n_lbd,n_int8+1:n_int8+n_lfd) = speye(n_lbd);

% intermediate summations
n_int9 = n_i+n_l+n_r+n_f+n_b+n_d+n_t+n_lf+n_rf+n_rb+n_lb+n_ld+n_rd+n_rt+n_lt+n_fd+n_bd+n_bt+n_ft+n_lfd+n_rfd+n_rbd+n_lbd;

T_per.s3(n_int9+1:n_int9+n_lft,n_int8+1:n_int8+n_lfd) = speye(n_lft);
T_per.s13(n_int9+n_lft+1:n_int9+n_lft+n_rft,n_int8+1:n_int8+n_lfd) = speye(n_rft);
T_per.s123(n_int9+n_lft+n_rft+1:n_int9+n_lft+n_rft+n_rbt,n_int8+1:n_int8+n_lfd) = speye(n_rbt);
T_per.s23(n_int9+n_lft+n_rft+n_rbt+1:n_int9+n_lft+n_rft+n_rbt+n_lbt,n_int8+1:n_int8+n_lfd) = speye(n_lbt);

% unsort rows to correspond to original DOF ordering
% (row "unsort" is important, but don't need to "unsort" columns...
%  I'm on the fence about keeping it in)
% T_per.s0   = T_per.s0(i_unsort1,i_unsort2);
% T_per.s1   = T_per.s1(i_unsort1,i_unsort2);
% T_per.s2   = T_per.s2(i_unsort1,i_unsort2);
% T_per.s3   = T_per.s3(i_unsort1,i_unsort2);
% T_per.s12  = T_per.s12(i_unsort1,i_unsort2);
% T_per.s13  = T_per.s13(i_unsort1,i_unsort2);
% T_per.s23  = T_per.s23(i_unsort1,i_unsort2);
% T_per.s123 = T_per.s123(i_unsort1,i_unsort2);


T_per.s0   = T_per.s0(i_unsort1,:);
T_per.s1   = T_per.s1(i_unsort1,:);
T_per.s2   = T_per.s2(i_unsort1,:);
T_per.s3   = T_per.s3(i_unsort1,:);
T_per.s12  = T_per.s12(i_unsort1,:);
T_per.s13  = T_per.s13(i_unsort1,:);
T_per.s23  = T_per.s23(i_unsort1,:);
T_per.s123 = T_per.s123(i_unsort1,:);