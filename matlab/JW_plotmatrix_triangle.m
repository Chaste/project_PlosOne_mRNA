function [h,ax,BigAx] = JW_plotmatrix_triangle(varargin)
%PLOTMATRIX Scatter plot matrix.
%   PLOTMATRIX(X,Y) scatter plots the columns of X against the columns
%   of Y.  If X is P-by-M and Y is P-by-N, PLOTMATRIX will produce a
%   N-by-M matrix of axes. PLOTMATRIX(Y) is the same as PLOTMATRIX(Y,Y)
%   except that the diagonal will be replaced by HIST(Y(:,i)).
%
%   PLOTMATRIX(...,'LineSpec') uses the given line specification in the
%   string 'LineSpec'; '.' is the default (see PLOT for possibilities).
%
%   PLOTMATRIX(AX,...) uses AX as the BigAx instead of GCA.
%
%   [H,AX,BigAx,P,PAx] = PLOTMATRIX(...) returns a matrix of handles
%   to the objects created in H, a matrix of handles to the individual
%   subaxes in AX, a handle to big (invisible) axes that frame the
%   subaxes in BigAx, a matrix of handles for the histogram plots in
%   P, and a matrix of handles for invisible axes that control the
%   histogram axes scales in PAx.  BigAx is left as the CurrentAxes so
%   that a subsequent TITLE, XLABEL, or YLABEL will be centered with
%   respect to the matrix of axes.
%
%   Example:
%       x = randn(50,3); y = x*[-1 2 1;2 0 1;1 -2 3;]';
%       plotmatrix(y)

%   Clay M. Thompson 10-3-94
%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.19.4.9 $  $Date: 2009/10/24 19:18:36 $

% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
%error(nargchk(1,3,nargs,'struct'));
nin = nargs;

sym = '.'; % Default scatter plot symbol.
dohist = 1;

if ischar(args{nin}),
    sym = args{nin};
    [l,c,m,msg] = colstyle(sym);
    if ~isempty(msg), error(msg); end
    nin = nin - 1;
end

% if nin==2, % plotmatrix(x,x) twice
    rows = size(args{2},2); cols = size(args{1},2);
    x = args{1}; y = args{1};
    w = args{2}; z = args{2};
    x_shared_indices = args{3}; y_shared_indices = args{3};
    w_shared_indices = args{4}; z_shared_indices = args{4};
% else
%     error('MATLAB:my_plotmatrix:InvalidLineSpec',...
%         'Invalid marker specification. Type ''help plot''.');
% end

if ndims(x)>2 || ndims(y)>2,
    error(id('InvalidXYMatrices'),'X and Y must be 2-D.')
end
if size(x,1)~=size(y,1) || size(x,3)~=size(y,3),
    error(id('XYSizeMismatch'),'X and Y must have the same number of rows and pages.');
end

% Create/find BigAx and make it invisible
BigAx = newplot(cax);
fig = ancestor(BigAx,'figure');
hold_state = ishold(BigAx);
set(BigAx,'Visible','off','color','none')

if any(sym=='.'),
    units = get(BigAx,'units');
    set(BigAx,'units','pixels');
    pos = get(BigAx,'Position');
    set(BigAx,'units',units);
    markersize = max(5,min(15,round(15*min(pos(3:4))/max(1,size(x,1))/max(rows,cols))));
else
    markersize = get(0,'DefaultLineMarkerSize');
end

% Create and plot into axes
ax = zeros(rows-1,cols-1);
pos = get(BigAx,'Position');
width = pos(3)/(cols-1);
height = pos(4)/(rows-1);
space = .02; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height];
m1 = size(y,1);
k1 = size(y,3);
m2 = size(z,1);
k2 = size(z,3);
xlim = -ones([rows-1 cols-1 2]);
ylim = -ones([rows-1 cols-1 2]);
BigAxHV = get(BigAx,'HandleVisibility');
BigAxParent = get(BigAx,'Parent');
paxes = findobj(fig,'Type','axes','tag','PlotMatrixScatterAx');
for i=rows-1:-1:1,
    for j=cols-1:-1:1,
    %for j=1:cols-1,
        if i >= j
            axPos = [pos(1)+(j-1)*width pos(2)+(rows-1-i)*height ...
                width*(1-space) height*(1-space)];
            findax = findaxpos(paxes, axPos);
            if isempty(findax),
                ax(i,j) = axes('Position',axPos,'HandleVisibility',BigAxHV,'parent',BigAxParent);
                set(ax(i,j),'visible','on');
            else
                ax(i,j) = findax(1);
            end
            
            %         pixels_across_subplot = 500;
            %         [plot_x plot_y plot_w plot_z plot_a plot_b] = make_new_set(...
            %             reshape(x(:,j,:),[m1 k1]), reshape(y(:,i,:),[m1 k1]),...
            %             reshape(w(:,j,:),[m2 k2]), reshape(z(:,i,:),[m2 k2]),...
            %             pixels_across_subplot);
            plot_x = reshape(x(:,j,:),[m1 k1]);
            plot_y = reshape(y(:,i+1,:),[m1 k1]);
            plot_w = reshape(w(:,j,:),[m2 k2]);
            plot_z = reshape(z(:,i+1,:),[m2 k2]);
            plot_a = [plot_x(x_shared_indices); plot_w(w_shared_indices)];
            plot_b = [plot_y(y_shared_indices); plot_z(z_shared_indices)];
            plot_x(x_shared_indices) = [];
            plot_y(y_shared_indices) = [];
            plot_w(w_shared_indices) = [];
            plot_z(z_shared_indices) = [];
            % [64 16 207] is the rgb "transparent" colour,
            % but not very striking
%             set(0,'DefaultAxesColorOrder',[0.65*[1 1 1]; ... % 'red'
%                 0.3*[1 1 1]; ... % blue
%                 0.1*[1 1 1]]); % yellow 229/256 255/256 0/256
            set(0,'DefaultAxesColorOrder',[[1 0 0]; ... % 'red'
                [0 0 1]; ... % blue
                [51/255 0 115/255]]); % yellow 229/256 255/256 0/256
            temp = plot(plot_x, plot_y,'.',...
                plot_w, plot_z,'.',...
                plot_a, plot_b,'.',...
                'parent',ax(i,j));
            
            hh(i,j,:) = temp(1:2)';
            
            set(temp,'markersize',markersize);
            set(ax(i,j),'xlimmode','auto','ylimmode','auto','xgrid','off','ygrid','off')
            xlim(i,j,:) = get(ax(i,j),'xlim');
            ylim(i,j,:) = get(ax(i,j),'ylim');
        end
    end
end

for i=rows-1:-1:1
    for j=cols-1:-1:1
        if i<j
            xlim(i,j,:) = get(ax(rows-1,j),'xlim');
            ylim(i,j,:) = get(ax(i, 1), 'ylim');
        end
    end
end

xlimmin = min(xlim(:,:,1),[],1);
xlimmax = max(xlim(:,:,2),[],1);
ylimmin = min(ylim(:,:,1),[],2);
ylimmax = max(ylim(:,:,2),[],2);

% Try to be smart about axes limits and labels.  Set all the limits of a
% row or column to be the same and inset the tick marks by 10 percent.
inset = .05;
dy = zeros(1,rows-1);
dx = zeros(1,cols-1);

% Go along diagonal
for i=1:rows-1
        set(ax(i,i),'ylim',[ylimmin(i,1) ylimmax(i,1)])
        set(ax(i,i),'xlim',[xlimmin(1,i) xlimmax(1,i)])
        dy(i) = diff(get(ax(i,i),'ylim'))*inset;
        dx(i) = diff(get(ax(i,i),'xlim'))*inset;
end

for i=1:rows-1,
    for j=1:cols-1,       
        if i >= j
            set(ax(i,j),'ylim',[ylimmin(i,1)-dy(i) ylimmax(i,1)+dy(i)])
            set(ax(i,j),'xlim',[xlimmin(1,j)-dx(j) xlimmax(1,j)+dx(j)])
        end
    end
end


for i=rows-1:-1:1,
    for j=cols-1:-1:1,
        if i >= j
            if j>=2
            set(ax(i,j),'yticklabel','');
            end
            if i<5
            set(ax(i,j),'xticklabel','');
            end
        end
    end
end
set(BigAx,'XTick',get(ax(rows-1,1),'xtick'),'YTick',get(ax(rows-1,1),'ytick'), ...
    'userdata',ax,'tag','PlotMatrixBigAx')
set(ax,'tag','PlotMatrixScatterAx');

% if dohist, % Put a histogram on the diagonal for plotmatrix(y) case
%     paxes = findobj(fig,'Type','axes','tag','PlotMatrixHistAx');
%     pax = zeros(1, rows);
%     for i=rows:-1:1,
%         axPos = get(ax(i,i),'Position');
%         findax = findaxpos(paxes, axPos);
%         if isempty(findax),
%             histax = axes('Position',axPos,'HandleVisibility',BigAxHV,'parent',BigAxParent);
%             set(histax,'visible','on');
%         else
%             histax = findax(1);
%         end
%         nbins = 100;
%         [nn,xx] = hist(reshape(y(:,i,:),[m1 k1]),nbins);
%         patches(i,:) = bar(histax,xx,nn,'hist');
%         %set(patches(i,:),'FaceColor','r','EdgeColor','none','facealpha',0.75);
%         set(patches(i,:),'FaceColor',0.55*[1 1 1],'EdgeColor','none','facealpha',1.0, 'edgealpha', 0.5);
%         hold on
%         [mm,zz] = hist(reshape(z(:,i,:),[m2 k2]),nbins);        
%         patches2=bar(histax,zz,mm,'hist');
%         %set(patches2,'FaceColor','b','EdgeColor','none','facealpha',0.75);
%         set(patches2,'facecolor',0*[1 1 1],'facealpha',0.75, 'EdgeColor', 'none', 'edgealpha', 0.5);
%         set(histax,'xtick',[],'ytick',[],'xgrid','off','ygrid','off');
%         set(histax,'xlim',[xlimmin(1,i)-dx(i) xlimmax(1,i)+dx(i)])
%         set(histax,'tag','PlotMatrixHistAx');
%         pax(i) = histax;  % ax handles for histograms
%     end
%     patches = patches';
% end

% Make BigAx the CurrentAxes
set(fig,'CurrentAx',BigAx)
if ~hold_state,
    set(fig,'NextPlot','replace')
end

% Also set Title and X/YLabel visibility to on and strings to empty
set([get(BigAx,'Title'); get(BigAx,'XLabel'); get(BigAx,'YLabel')], ...
    'String','','Visible','on')

if nargout~=0,
    h = hh;
end

function str=id(str)
str = ['MATLAB:my_plotmatrix:' str];

function findax = findaxpos(ax, axpos)
tol = eps;
findax = [];
for i = 1:length(ax)
    axipos = get(ax(i),'Position');
    diffpos = axipos - axpos;
    if (max(max(abs(diffpos))) < tol)
        findax = ax(i);
        break;
    end
end

