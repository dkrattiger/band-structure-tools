function V = Si_Pseudo(r,a0);

Ry2Jo = 2.1799e-18;
Ry2Jo = 1;

% Pseudopotential Form Factors
V3S = -0.21*Ry2Jo;
V8S = 0.04*Ry2Jo;
V11S = 0.08*Ry2Jo;

% lattice vectors
% a0 = 0.5431e-9;
% a0 = 10.261;
% a0 = 1;
% a0 = 5.43;

% a1 = [0;1;1];
% a2 = [1;0;1];
% a3 = [1;1;0];

% reciprocal lattice vectors
b1 = [-1;1;1];
b2 = [1;-1;1];
b3 = [1;1;-1];

% atom locations
tau = a0*[1/8;1/8;1/8];
r = r-tau*ones(1,size(r,2));

% Vectorized Pseuodpotential calculation
ijk_lim = 2; % this limit encompasses all points such that |G|<=11
ijk_vec = [-ijk_lim:ijk_lim];
[ii,jj,kk] = ndgrid(ijk_vec,ijk_vec,ijk_vec);

ii = ii(:);
jj = jj(:);
kk = kk(:);

G = [b1,b2,b3]*[ii,jj,kk]';
Gnorm = sum(G.*G);

Valt=0;

i_G3 = Gnorm==3;
i_G8 = Gnorm==8;
i_G11 = Gnorm==11;

Gtau = G'*tau;

% size(G)
% size(r)
% i_G3

Gr = G'*r;

% size(V3S*cos((2*pi/a0)*Gtau(i_G3))')
% size(exp(-1i*(2*pi/a0)*Gr(i_G3,:)))

Valt = Valt + V3S*cos((2*pi/a0)*Gtau(i_G3))'*exp(-1i*(2*pi/a0)*Gr(i_G3,:));
Valt = Valt + V8S*cos((2*pi/a0)*Gtau(i_G8))'*exp(-1i*(2*pi/a0)*Gr(i_G8,:));
Valt = Valt + V11S*cos((2*pi/a0)*Gtau(i_G11))'*exp(-1i*(2*pi/a0)*Gr(i_G11,:));

% taking just the real part seems to be fine (results didn't change)
Valt = real(Valt);
V = Valt;