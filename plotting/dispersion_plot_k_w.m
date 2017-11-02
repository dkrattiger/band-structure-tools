function [h] = dispersion_plot_k_w(omegas,kappas,options)

% set line width parameter
linewidths = 1;

% set tolerance for imaginary/real cutoffs
tol = 1e-2;

%% set default plot settings
% ======================================================================= %
dfcolor = get(groot,'factoryAxesColorOrder');% matlab's builtin color order

defaults.Colors             = mat2cell(dfcolor,ones(1,size(dfcolor,1)),[3]); 
defaults.LineStyles         = {'-','--',':','-.'};
defaults.LineWidths         = {1};
defaults.MarkerEdgeColors   = {'auto'};  
defaults.MarkerFaceColors   = {'none'};   
defaults.Markers            = {'none'};
defaults.SeparateSubplots   = true;
defaults.splitImag          = false;
defaults.ThreeD             = false;
% defaults.BZlabels           = {};   % strings for x-tick-labels
% defaults.BZLabelLocations   = [];   % Where to place labels (and dashed lines)
defaults.legendstrings      = {};   % cell array of strings containing legend names
defaults.add_labels         = true;

% start with default settings and overwrite any that are given in "options"
if nargin>=3
    options = setstructfields(defaults,options);
else
    options = defaults;
end

%% Check the arguments "kappas" and "omega"
% ======================================================================= %

% if kappas is not a cell array, make it one
if ~iscell(kappas)
    kappas = {kappas};
end


% if omega is not a cell array, make it one
if ~iscell(omegas)
    omegas = {omegas};
end

% if just a single vector is provided for omega, copy it for each instance
% in kappas
if length(omegas) == 1
    for i = 2:length(kappas)
        omegas{i} = omegas{1};
    end
end



%% loop through each set of results being plotted
% ======================================================================= %
for i = 1:length(kappas)
    
    % extract current set of kappa and omega solutions from cell array
    kappa = kappas{i};
    omega = omegas{i};
    
    % find maximum kappa solution
    max_kap = max(max(abs(real(kappa))));
    
    % create logical index of real kappas
    i_real = abs(real(kappa))>max_kap*tol | ...
        (abs(real(kappa))<max_kap*tol & abs(imag(kappa))<max_kap*tol);
    
    % create logical index of imaginary kappas
    i_imag = abs(imag(kappa))>max_kap*tol;
    
    % create logical index of kappas at edge of BZ
    i_X = abs(max_kap-abs(real(kappa)))<max_kap*tol & i_imag;
        
    % create array to plot real dispersion
    kappa_plot{2} = abs(real(kappa));
    kappa_plot{2}(i_imag) = nan;
    
    % create array to plot imaginary dispersion
    kappa_plot{1} = -abs(imag(kappa));
    kappa_plot{1}(i_real) = nan;
    
    % create array to plot imaginary dispersion at BZ edge
    kappa_plot{3} = abs(imag(kappa));
    kappa_plot{3}(~i_X) = nan;
    
    % plot and format dispersion curves
    %SeparateSubplots = false;
    if options.SeparateSubplots
            
        for j = 1:3

            subplot(1,3,j)
            h{i,j} = plot((kappa_plot{j}),omega,'linewidth',linewidths);hold on

            if j == 2
                xlim([0,max_kap])
            end
            ylim(omega([1,end]))

            % add labels
            if options.add_labels
                xsublabels = {{'imaginary';'wavevectors'},...
                              {'real';'wavevectors'},...
                              {'complex';'wavevectors'}};
                xlabel(xsublabels{j});
                if j==1
                    ylabel('Frequency, \Omega');
                else
                    set(gca,'yticklabel',[])
                end
            end
        end
        
        
        
    else % single plot axes
        
        kappa_plot{3} = kappa_plot{3}+max_kap;
     
        % loop through imaginary, real, complex
        for j = 1:3
            h{i,j} = plot((kappa_plot{j}),omega,'linewidth',linewidths);hold on
        end
        
        % set axes limits
        ylim(omega([1,end]))        

        % plot lines showing edge of real/imaginary wavenumber regions
        set(gca,'xtick',[0,max_kap])
        set(gca,'xgrid','on')           
        
    end
    
    for j = 1:3
                    
        % set line color
        i_color = rem(i-1,size(options.Colors,1))+1;
        set(h{i,j},'color', options.Colors{i_color});
        
        % set line style
        i_linestyle = rem(i-1,length(options.LineStyles))+1;
        set(h{i,j},'linestyle', options.LineStyles{i_linestyle});
        
        % set linewidths
        i_linewidth = rem(i-1,length(options.LineWidths))+1;
        set(h{i,j},'linewidth', options.LineWidths{i_linewidth});
        
        % set marker type
        i_marker = rem(i-1,length(options.Markers))+1;
        set(h{i,j},'marker', options.Markers{i_marker});
        
        % set marker edge color
        i_markerEdge = rem(i-1,length(options.MarkerEdgeColors))+1;
        set(h{i,j},'markeredgecolor', options.MarkerEdgeColors{i_markerEdge});
        
        % set marker edge color
        i_markerFace = rem(i-1,length(options.MarkerFaceColors))+1;
        set(h{i,j},'markerfacecolor', options.MarkerFaceColors{i_markerFace});
    end
    
    % define legend handle vector
    legendvec(i) = h{i,1}(1);
end

% add legend
if ~isempty(options.legendstrings)
    legend(legendvec,options.legendstrings,'location','northeast');
end
hold off