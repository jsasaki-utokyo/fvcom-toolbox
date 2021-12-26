function plot_TIDE_current(Num, strds, lon, lat, ua, va, DateString)

% Import working dic
if ismac
    addpath('/Users/yulong/Documents/GitHub/fvcom-toolbox/m_map/');
    addpath('/Users/yulong/Documents/GitHub/fvcom-toolbox/t_tide/');
    addpath('/Users/yulong/Documents/GitHub/fvcom-toolbox/tidal_ellipse/');
elseif isunix
    addpath('/home/usr0/n70110d/github/fvcom-toolbox/m_map/');
    addpath('/home/usr0/n70110d/github/fvcom-toolbox/t_tide/');
    addpath('/home/usr0/n70110d/github/fvcom-toolbox/tidal_ellipse/');
elseif ispc
    addpath('C:/Users/Yulong WANG/Documents/GitHub/fvcom-toolbox/m_map/');
    addpath('C:/Users/Yulong WANG/Documents/GitHub/fvcom-toolbox/t_tide/');
    addpath('C:/Users/Yulong WANG/Documents/GitHub/fvcom-toolbox/tidal_ellipse/');
end

% test code
% load('tidal_ellipse_ncfile.mat');
% load('tidal_ellipse_out2.mat');
% Num = 20;
% strds = 500;
% lon = ncfile.lonc;
% lat = ncfile.latc;
% ua = ncfile.ua(:,1);
% va = ncfile.va(:,1);
% time = ncfile.time(1);

% DateString = datestr(double(time), 31);

% Create figure
figure1 = figure(Num);
clf

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
m_proj('transverse mercator','lon',[139.1,140.4],'lat',[34.65,35.75]);
for i = [1:1:2000,2001:1:3000,3001:1:3500,3401:5:length(lon)]
    [llx,lly] = m_ll2xy(lon(i),lat(i),'clip','off');
    if va(i) > 0
        quiver(llx,lly,ua(i)/strds,va(i)/strds,'r-');
    elseif va(i) < 0
        quiver(llx,lly,ua(i)/strds,va(i)/strds,'b-');
    end
    hold on;
end
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','k'); 
m_grid('linest','none','box','on','xtick',9);
% Create ylabel
ylabel('Latitude (degrees)');
% Create xlabel
xlabel('Longtitude (degrees)');
% Create title
title(['Tidal current: ', DateString]);
set(gcf,'position',[20*10,20*10,600,600])

