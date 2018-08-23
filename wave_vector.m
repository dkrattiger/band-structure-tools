function [kappa,kappa_plot] = wave_vector(symPts,refLevel,R)
%
% Dimitri Krattiger
%
% Description
% ===========
% This function produces a discrete set of wave vectors that trace the edge
% of the Brillouin Zone
%
% inputs
% ======
% sym_pts   = cell array of strings specifying the BZ symmetry points
%
% refLevel  = refinement level incrementing this by 1 roughly doubles the 
%             number of k-pts
%
% R         = lattice vectors
%
% outputs
% =======
% kappa      = discrete set of wave vectors
% 
% kappa_plot = gives spacing between k-points in kappa (good for plotting
%              band structure diagram)


% number of dimensions and periodicity directions
[n_d,n_p] = size(R);

% compute reciprocal lattice vectors (B) using matrix representation
B = 2*pi*inv(R(1:n_p,1:n_p))';

if n_p == 1
    
    % define high symmetry points
    sym_pt_struct.Gamma = B*0;
    sym_pt_struct.X = B*1/2;
    
elseif n_p == 2
    
    % define high symmetry points
    sym_pt_struct.Gamma = B*[0;0];
    sym_pt_struct.X = B*[1/2;0];
    sym_pt_struct.Y = B*[0;1/2];
    sym_pt_struct.M = B*[1/2;1/2];
    
elseif n_p == 3
    
    % define high symmetry points
    sym_pt_struct.Gamma = B*[0;0;0];
    sym_pt_struct.X = B*[1/2;0;0];
    sym_pt_struct.M = B*[1/2;1/2;0];
    sym_pt_struct.R = B*[1/2;1/2;1/2];
    sym_pt_struct.M_10 = (sym_pt_struct.M)/10;
    sym_pt_struct.X_5 = (sym_pt_struct.X)/5;
end
    
% collect BZ symmetry points into an array
n_sym_pts = length(symPts);



for i = 1:n_sym_pts
    if symPts{i}(1) == '\';
        symPts{i}(1) = [];
    end
end

for i = 1:n_sym_pts
    sym_pt_array(:,i) = sym_pt_struct.(symPts{i});
end

% number of k-points
n_kap = (n_sym_pts-1)*(2^refLevel)+1;

% number of k-points per BZ segment
n_kap_seg = (n_kap-1)/(n_sym_pts-1)+1;

% define wave-vector for every k-point
kappa = zeros(n_d,n_kap);

for i = 1:n_sym_pts-1
    
    % index for current wavevector segment   
    i_kap = ((i-1)*(n_kap_seg-1)+1):(i*(n_kap_seg-1)+1);
    
    kappa(1:n_p,i_kap) = ...
        sym_pt_array(:,i)*ones(1,n_kap_seg)+(sym_pt_array(:,i+1)-sym_pt_array(:,i))*...
        linspace(0,1,n_kap_seg);
end

kappa_plot = zeros(1,n_kap);
for i = 2:n_kap
    kappa_plot(i) = kappa_plot(i-1)+norm(kappa(:,i)-kappa(:,i-1));
end
kappa_plot = kappa_plot;