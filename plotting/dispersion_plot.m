function [h_disp] = dispersion_plot(kappa,omegas,varargin)

if ~iscell(omegas)
    omegas = {omegas};
end

% check hold status
hold_on = ishold;

% % % Get legend info for current plot
h_legend = findobj(gca,'type','legend');
old_legendstrings = get(h_legend,'string');
old_legendvec = get(h_legend,'userdata');

% track how many sets of dispersion curves have been plotted and update
% this in user-data field of figure object
n_om_sets_old = get(gca,'userdata');
if isempty(n_om_sets_old)
    n_om_sets_old = 0;
end
n_om_sets = n_om_sets_old + length(omegas);
set(gca,'userdata',n_om_sets);

% set line width for all curves
linewidths = 2;

% number of k-points
kappa = kappa(:);
n_kap = length(kappa);

% get default color and linstyle orders
colors = get(gca,'ColorOrder');
linestyles = {'-',':','-.','-',':','-.','-',':','-.'};
% linestyles = {'-'};


% loop through each data set and plot and format dispersion curves
for i = 1:length(omegas)
    h_disp{i} = plot(kappa,omegas{i},'linewidth',linewidths);hold on
    set(h_disp{i},'color', colors(rem(i+n_om_sets_old-1,size(colors,1))+1,:));
    set(h_disp{i},'linestyle', linestyles{rem(i+n_om_sets_old-1,length(linestyles))+1});
    legendvec(i) = h_disp{i}(1);
end

%set xlimit
xlim([kappa(1),kappa(n_kap)]);

% plot and label high symmetry points
if nargin>=3
    sym_pts = varargin{1};
    if ~isempty(sym_pts)
        n_sym_pts = length(sym_pts);

        % set and label ticks
        i_tick = linspace(1,n_kap,n_sym_pts);
        set(gca,'xtick',kappa(i_tick))   
        set(gca,'xticklabel',sym_pts)

        % add dashed lines at sym pts
        ylimits = ylim;
        hold on
        plot([1;1]*kappa(i_tick(2:n_sym_pts-1)).',ylimits'*ones(1,n_sym_pts-2),'k:')
        hold off
    end
end

% return to hold status that we started with
if hold_on
    hold on
else
    hold off
end

% add legend
if nargin>=4
    legendstrings = [old_legendstrings,varargin{2}];
    legendvec = [old_legendvec,legendvec];
    if ~isempty(legendstrings)
        h_legend = legend(legendvec,legendstrings,'location','northeast');
        set(h_legend,'userdata',legendvec);
    end
end
xlabel('Wave vector');ylabel('Frequency, \Omega');%title('Dispersion')

% change colors
if nargin>=5
    colors = varargin{3};
    for i = 1:length(omegas)
        set(h_disp{i},'color',colors{i})
    end
end

% change linestyles
if nargin>=6
    linestyles = varargin{4};
    for i = 1:length(omegas)
        set(h_disp{i},'linestyle',linestyles{i})
    end
end

% change markers
if nargin>=7
    linestyles = varargin{5};
    for i = 1:length(omegas)
        set(h_disp{i},'marker',markers{i})
    end
end