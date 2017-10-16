function [PHI,L,iWindow] = eigsFull(K,M,varargin)

%% default options
% ======================================================================= %
defaults.verbose = true;
defaults.k = 30;
defaults.plot = false;
%defaults.p_eig = 3*k; (can't be defined before options are resolved)

if nargin>=3
    options = varargin{1};
else
    options.dummy = [];
end

options = setstructfields(defaults,options);

% Problem size
n = size(K,1);

% number of eigenvalues to compute per window
k = options.k; %number of

% specify number of Lanczos basis vectors to be used in eigs. The default
% used by eigs is p=2k (p=2k+1 for hermitian). I observed that this caused the
% last few computed eigenvalues to be inaccurate when running a Silicon
% atomic structure unit cell calculation. It works well for p=3k, but I
% haven't investigated this further.
options.p_eig = min(3*k,n);


% create new options structue to be used by eigs() function
options_eigs.p = options.p_eig;

%% compute low and high frequency eigenvalues first
% ======================================================================= %
[PHIcell{1},Ltemp] = eigs(K,M,k,'sm',options_eigs);
Lcell{1} = diag(Ltemp);
Lmaxs(1) = max(Lcell{1}); Lmins(1) = min(Lcell{1});

[PHIcell{2},Ltemp] = eigs(K,M,k,'lm',options_eigs);
Lcell{2} = diag(Ltemp);
Lmaxs(2) = max(Lcell{2}); Lmins(2) = min(Lcell{2});

% compute eigenalue range
eigRange = (max(Lmaxs)-min(Lmins));

% determine sorting index for solution cells
Lcenters = (Lmaxs+Lmins)/2;
[~,isortcell] = sort(Lcenters);

% determine gaps in eigenvalue coverage
gaps = Lmins(isortcell(2:end))-Lmaxs(isortcell(1:end-1));
[gapmax,imax] = max(gaps);


%% compute additional eigenvalue cell solutions
% ======================================================================= %
% if the maximum gap size is greater than or equal to zero, use the center
% of the gap as the next centering frequency

% tolerance to check degeneracy
degenTol = eigRange*1e-9; % try this tolerance for now
% degenTol = eigRange*1e-8; % try this tolerance for now


% plot eigenvalue calculation and windows
if options.plot
    figure(123);clf;grid on; hold on
    for j = 1:2
        plot(Lcell{j}/Lmaxs(2),0*ones(1,k),'k.')
    end
    plot([Lmins;Lmaxs]/Lmaxs(2),0*ones(2,2),'-','color',[0,0,1,0.5],'linewidth',10)
end

cellcount = 2;

% iterate until there is a minimum amount of overlap between all windows
while gapmax >= -degenTol
    
    tic
    % increment cell counter
    cellcount = cellcount+1;
    
    % compute range of current gap and center of current gap
    Lrange = Lmaxs(isortcell(imax))-Lmins(isortcell(imax+1));
    Lcenter = Lmins(isortcell(imax+1)) + Lrange*(1/2);
    
    % compute eigenvalues centered in current gap
    try
        [PHIcell{cellcount},Ltemp] = eigs(K,M,k,Lcenter,options_eigs);
        Lcell{cellcount} = diag(Ltemp);
        Lmaxs(cellcount) = max(Lcell{cellcount}); 
        Lmins(cellcount) = min(Lcell{cellcount});
    catch        
        % if centering frequency happens to fall on an eigenvalue, perturb
        % it by a small amount (average distance btwn eigvals)
        Lcenter = Lcenter + eigrange/n;
        [PHIcell{cellcount},Ltemp] = eigs(K,M,k,Lcenter,options_eigs);
        Lcell{cellcount} = diag(Ltemp);
        Lmaxs(cellcount) = max(Lcell{cellcount}); 
        Lmins(cellcount) = min(Lcell{cellcount});
    end
    
    % recompute sorting index for solution cells
    Lcenters = (Lmaxs+Lmins)/2;
    [~,isortcell] = sort(Lcenters);

    % determine gaps in eigenvalue coverage
    gaps = Lmins(isortcell(2:end))-Lmaxs(isortcell(1:end-1));
    [gapmax,imax] = max(gaps);
    
    % plot eigenvalue calculation and windows        
    if options.plot
        for j = 1:cellcount
            plot(Lcell{j}/Lmaxs(2),-(cellcount-2)*ones(1,k),'k.')
        end
        plot([Lmins;Lmaxs]/Lmaxs(2),-(cellcount-2)*ones(2,cellcount),'-','color',[0,0,1,0.5],'linewidth',10)
        plot(Lcenter/Lmaxs(2),-(cellcount-3),'o','color','k','markersize',12,'markerfacecolor','r')
    end
    
    if options.verbose
        fprintf('\t%i cells computed, total eigenvalue range = %4.2f, max gap size = %4.2e, loop time: %4.2f\n',...
            cellcount,eigRange,gapmax,toc);
    end
end
drawnow

%% Sort through computed solutions and discard duplicates
% ======================================================================= %

L = nan(n,1);
PHI = nan(n,n);
iWindow = zeros(n,1);

% start off by keeping all solutions from the first cell
[L(1:k),isorteig] = sort(Lcell{isortcell(1)});
PHI(:,1:k) = PHIcell{isortcell(1)}(:,isorteig);
iWindow(1:k) = 1;

eigcount = k;
for i = 2:cellcount
    
    j = isortcell(i);
    [L2,isorteig] = sort(Lcell{j});
    PHI2 = PHIcell{j}(:,isorteig);
    
    % last added value of L (should be largest)
    Lmax = L(eigcount);
    
    % compute number of terms in new eigenvalue window that overlap with
    % old set
    
    % maximum number of terms that overlap
    n_ovlp_max = sum(L2<(Lmax+2*degenTol));
    
    % check all possible overlaps to find the one that gives smallest
    % difference
    diffnorm = zeros(1,n_ovlp_max-1);
    for j = 1:n_ovlp_max
        
        L_ovlp = L((eigcount-j+1):eigcount);
        L2_ovlp = L2(1:j);
        
        diffnorm(j) = sum(abs(L_ovlp-L2_ovlp))/j;
    end
    
    [~,imin] = min(diffnorm);
    n_ovlp = imin;
    
   
    if n_ovlp_max==1
        warning(['overlap of just 1 occurred between two eigenvalue ',... 
        'sets. This Indicates that the error in computed solutions ',...
        'between two sets may be similar to the size of the overlap ',...
        'tolerance.'])
    elseif n_ovlp_max==0
        warning(['no overlap detected... no idea why this would happen.',...
            'Manual fix applied, but do not trust this solution'])
        n_ovlp = 0;
    end
    
    % add eigenvalues in L2 that don't overlap with L into L
    i_keep2 = (n_ovlp+1):k;
    L(eigcount + (1:length(i_keep2))) = L2(i_keep2);
    PHI(:,eigcount + (1:length(i_keep2))) = PHI2(:,i_keep2);
    iWindow(eigcount + (1:length(i_keep2))) = i;
    
    % increment eigenvalue counter
    eigcount = eigcount+length(i_keep2);
        
    if options.verbose
        fprintf('%i overlapping eigenvalues\n',n_ovlp)
    end
end
