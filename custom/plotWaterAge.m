function plotWaterAge(ncfile,clrbar)
% plotWaterAge("test_0001.nc",0.1)
%     clear all;
%     clc;

    if ismac    % On Mac
        basedir = '/Users/yulong/GitHub/';
        addpath([basedir,'fvcomtoolbox/']);
        addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
        addpath([basedir,'fvcomtoolbox/utilities/']);
        addpath([basedir,'fvcomtoolbox/custom/']);
    elseif isunix       % Unix?
        basedir = '/home/usr0/n70110d/github/';
        addpath([basedir,'fvcomtoolbox/']);
        addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
        addpath([basedir,'fvcomtoolbox/utilities/']);
        addpath([basedir,'fvcomtoolbox/custom/']);
    elseif ispc     % Or Windows?
        basedir = 'C:/Users/Yulong WANG/Documents/GitHub/';      
        % basedir = 'C:/Users/Yulong/Documents/GitHub/';      
        addpath([basedir,'fvcomtoolbox/']);
        addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
        addpath([basedir,'fvcomtoolbox/utilities/']);
        addpath([basedir,'fvcomtoolbox/custom/']);
    end

    % ncfile = "test_0001.nc";
    % clrbar = 0.1;

    nc.dye       = ncread(ncfile,'DYE');
    nc.dye_age   = ncread(ncfile,'DYE_AGE');
    nc.time      = ncread(ncfile,'time');
    nc.Times     = ncread(ncfile,'Times');
    nc.Times     = nc.Times';

    nc.nv        = ncread(ncfile,'nv');
    nc.nv        = double(nc.nv);

    nc.lon       = ncread(ncfile,'lon');
    nc.lon       = double(nc.lon);
    nc.lat       = ncread(ncfile,'lat');
    nc.lat       = double(nc.lat);
    nc.lonlat    = [nc.lon, nc.lat];
    nc.lonlat    = double(nc.lonlat);
    nc.tri_lonlat= triangulation(nc.nv,nc.lonlat);
    nc.tri       = nc.tri_lonlat.ConnectivityList;
    nc.tri(:,[2 3]) =nc.tri(:,[3 2]);

    water_mask = squeeze(nc.dye(:,1,:)<=1E-9);
    water_age = squeeze(nc.dye_age(:,1,:)./nc.dye(:,1,:));
    water_age(water_mask)=nan;
    water_age(water_age<0)=nan;
    water_age_day = water_age/24;

    fig = figure(01);
    for tt = 1:length(nc.time)
        plotMesh(01, [nc.lon, nc.lat], nc.tri, water_age_day(:,tt),...
            'time',nc.Times(tt,1:19),...
            'dye age',{'Water age (day)',[0,clrbar]});
        drawnow
        frame = getframe(fig);
        im{tt} = frame2im(frame);
    end
    filename = '03_water_age.gif'; 
    for tt = 1:length(nc.time)
        [A,map] = rgb2ind(im{tt},256);
        if tt == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
        end
    end