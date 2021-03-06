function plot_TS_surface(Num, XData1, YData1, CData1, X1, Y1, S1, C1, X2, Y2, X3, Y3, Name)
%CREATEFIGURE(ZData1, YData1, XData1, CData1, X1, Y1, S1, C1, X2, Y2, X3, Y3)
%  NUM: figure number
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  XDATA1:  surface xdata
%  CDATA1:  surface cdata
%  X1:  scatter x
%  Y1:  scatter y
%  S1:  scatter s
%  C1:  scatter c
%  X2:  vector of x data
%  Y2:  vector of y data
%  X3:  vector of x data
%  Y3:  vector of y data

%  Auto-generated by MATLAB on 23-May-2019 15:13:12

% Create figure
figure1 = figure(Num);
clf

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create surface
% surface('Parent',axes1,'ZData',ZData1,'YData',YData1,'XData',XData1,...
%     'AlignVertexCenters','on',...
%     'EdgeColor','none',...
%     'CData',CData1);
pcolor(XData1, YData1, CData1);
shading flat

% Create scatter
scatter(X1,Y1,S1,C1,'MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0]);

% Create plot
plot(X2,Y2,'Marker','square','LineStyle','none','Color',[1 0 0]);

% Create plot
plot(X3,Y3,'Marker','o','LineStyle','none','Color',[1 1 1]);

% Create ylabel
ylabel('Latitude (degree)');

% Create xlabel
xlabel('Longtitude (degree)');

% Create title
if strcmpi(Name, 'temperature')
    item = 'temperature (^{\circ}C)';
else
    item = 'salinity (PSU)';
end
title(['The 1st layer HYCOM ', item]);

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[139.1 140.4]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[34.65 35.75]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'CLim',[min(min(min(CData1))) max(max(max(CData1)))],'DataAspectRatio',...
    [1.2 1 1]);
set(gcf,'position',[Num*10,Num*10,800,600])
% Create colorbar
colorbar('peer',axes1);

