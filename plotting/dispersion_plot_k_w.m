function [h] = dispersion_plot_k_w(omega,kappas,options)

% set line width parameter
linewidths = 1;

% set tolerance for imaginary/real cutoffs
tol = 1e-2;

%% set default plot settings
% ======================================================================= %
dfcolor = get(groot,'factoryAxesColorOrder');% matlab's builtin color order

defaults.ColorOrder         = mat2cell(dfcolor,ones(1,size(dfcolor,1)),[3]); 
defaults.LineStyleOrder     = {'-','--',':','-.'};
defaults.LineWidthOrder     = {1};
defaults.MarkerEdgeColor    = {'auto'};  
defaults.MarkerFaceColor    = {'none'};   
defaults.MarkerOrder        = {'none'};
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

%% loop through each set of results being plotted
% ======================================================================= %
for i = 1:length(kappas)
    
    % extract current set of kappa solutions from cell array
    kappa = kappas{i};
    
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
        i_color = rem(i-1,size(options.ColorOrder,1))+1;
        set(h{i,j},'color', options.ColorOrder{i_color});
        
        % set line style
        i_linestyle = rem(i-1,length(options.LineStyleOrder))+1;
        set(h{i,j},'linestyle', options.LineStyleOrder{i_linestyle});
        
        % set linewidths
        i_linewidth = rem(i-1,length(options.LineWidthOrder))+1;
        set(h{i,j},'linewidth', options.LineWidthOrder{i_linewidth});
        
        % set marker type
        i_marker = rem(i-1,length(options.MarkerOrder))+1;
        set(h{i,j},'marker', options.MarkerOrder{i_marker});
        
        % set marker edge color
        i_markerEdge = rem(i-1,length(options.MarkerEdgeColor))+1;
        set(h{i,j},'markeredgecolor', options.MarkerEdgeColor{i_markerEdge});
        
        % set marker edge color
        i_markerFace = rem(i-1,length(options.MarkerFaceColor))+1;
        set(h{i,j},'markerfacecolor', options.MarkerFaceColor{i_markerFace});
    end
    
    % define legend handle vector
    legendvec(i) = h{i,1}(1);
end

% add legend
if ~isempty(options.legendstrings)
    legend(legendvec,options.legendstrings,'location','northeast');
end
hold off