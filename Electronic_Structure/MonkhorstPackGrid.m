function [k_ibz,w_ibz,varargout] = MonkhorstPackGrid(b1,b2,b3,nk,varargin)
% function [kpt,i_ibz,w_ibz] = MonkhorstPackGrid(b1,b2,b3,nk,varargin)
if nargin>=5
    options = varargin{1};
else
    options.dummy = [];
end

%% Check what inputs are given and set rest to default
% ======================================================================= %

% set_default values for options
defaults.grid_type      = 'GammaCentered'; % Monkhorst Pack
defaults.resort         = 'ks'; % Monkhorst Pack

% if ~exist('options')
%     options.dummy = 0;
% end
options = setstructfields(defaults,options);

%% 
% expand nk to 3x1 vector if necessary
if length(nk) == 1
    nk = nk*[1,1,1];
end

% create (nk x nk x nk) grid of k-points
switch options.grid_type
    case 'MonkhorstPack'
        
        % MP grid shift
        s = [0,0,0];
        use_shift = true;
        if use_shift
            s(rem(nk,2) == 0) = 1./(2*nk(rem(nk,2) == 0));
        end
        
        [i,j,k] = ndgrid(1:nk(1),1:nk(2),1:nk(3));
        i = i(:)'; j = j(:)'; k = k(:)';
%         kpt1 = [b1,b2,b3]*(([i;j;k]-nk(1)/2-1/2)/(nk)-s);
        kpt = [b1,b2,b3]*[(i-nk(1)/2-1/2)/nk(1)-s(1);...
                          (j-nk(2)/2-1/2)/nk(2)-s(2);...
                          (k-nk(3)/2-1/2)/nk(3)-s(3)];
       
    case 'GammaCentered'

        [i,j,k] = meshgrid(0:(nk(1)-1),0:(nk(2)-1),0:(nk(3)-1));
        i = i(:)'; j = j(:)'; k = k(:)';   
        kpt = [b1,b2,b3]*[i/nk(1);j/nk(2);k/nk(3)];
end

%% Sort k-points in specific way to get consistency across grid types

switch options.resort
    case 'bs'
        % create sorting array to define which points will be selected
        % 1st column (overall magnitude minimized)
        % 2nd column (b1 component of k maximized)
        % 3rd column (b2 component of k maximized)
        % 2nd column (b3 component of k maximized)
        ijk = [b1,b2,b3]\kpt;
        sorting_array = [sqrt(sum(kpt.^2))', -ijk(1,:)', -ijk(2,:)', -ijk(3,:)'];
        roundplace = 3;
        [~,isort] = sortrows(round(sorting_array*10^roundplace)/(10^roundplace));
    case 'ks'
        % create sorting array to define which points will be selected
        % 1st column (overall magnitude minimized)
        % 2nd column (kx maximized)
        % 3rd column (ky maximized)
        % 2nd column (kz maximized)
        sorting_array = [sqrt(sum(kpt.^2))', -kpt(1,:)', -kpt(2,:)', -kpt(3,:)'];
        roundplace = 3;
        [~,isort] = sortrows(round(sorting_array*10^roundplace)/(10^roundplace));
    case 'none'
        isort = 1:prod(nk);
end

% reorder kpts
kpt = kpt(:,isort);

%% Perform sym. operations and compare resulting grid with original to find 
%% symmetric points

tol = max(abs([b1;b2;b3]))*1e-3;
Lib_array = zeros(prod(nk),48);
for reflx = 0:1
    for refly = 0:1
        for reflz = 0:1

            p = sortrows(perms(1:3));
            np = size(p,1);
            for i = 1:size(p,1)
                
                kpt2 = kpt;
                
                if reflx
                    kpt2(1,:) = -kpt2(1,:);
                end
                if refly
                    kpt2(2,:) = -kpt2(2,:);
                end
                if reflz
                    kpt2(3,:) = -kpt2(3,:);
                end
                kpt2 = kpt2(p(i,:),:);
                
                ijk2 = [b1,b2,b3]\kpt2;
                tol_shift = 1e-3;
                
                % Map each point back to original sampled volume 
                switch options.grid_type                   
                    case 'MonkhorstPack'
                        ijk2 = mod(ijk2+1/2+tol_shift,1)-1/2-tol_shift;
                        kpt2 = [b1,b2,b3]*ijk2;
                    case 'GammaCentered'
                        ijk2 = mod(ijk2+tol_shift,1)-tol_shift;
                        kpt2 = [b1,b2,b3]*ijk2;
                end
                
                
                [Lia,Lib] = ismember_tol(kpt2',kpt','rows','tolerance',tol);
                Lib_array(:,2*2*np*reflx + 2*np*refly + np*reflz + i)  = Lib;
                
            end
        end
    end
end

% set zero values in Lib array (no match) to be one higher than the number
% of k-points
Lib_array(Lib_array==0) = prod(nk)+1;

symmetries = unique(sort(Lib_array,2),'rows');
w = zeros(size(symmetries,1),1);
for i = 1:size(symmetries,1)
    
    i_rep{i} = unique(symmetries(i,symmetries(i,:)~=0));
    i_exp(i_rep{i}) = i;
    w(i) = length(i_rep{i});
    
end
w_ibz = w/prod(nk);
i_ibz = symmetries(:,1);

% reduced set of k points
k_ibz = kpt(:,i_ibz);

%% specify outputs
if nargout >= 3
    varargout{1} = kpt;
end
if nargout >= 4
    varargout{2} = i_ibz;
end
if nargout >= 5
    varargout{3} = i_exp;
end

function [Lia,Lib] = ismember_tol(A,B,varargin)

% ismember function that determines set membership but within some
% tolerance


% find whether a tolerance optional argument is given 
% (if so remove it from the optiona argument list and save it)
option_args = varargin;
for j = (nargin-3):-1:1
    if strcmpi(varargin{j},'tolerance')
        tol = option_args{j+1};
        option_args([j,j+1]) = [];
    end
end

% use tolerance argument to determine number of rounding digits to use
if exist('tol')
    rd = 10^(-ceil(log10(tol)));
else
    rd = 1;
end

[Lia,Lib] = ismember(round(A*rd)/rd,round(B*rd)/rd,option_args{:});
