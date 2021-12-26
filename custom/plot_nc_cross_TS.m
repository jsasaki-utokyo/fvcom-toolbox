function plot_nc_cross_TS(Num, Mobj, ncfile, timestep, crs_node, vartoplot)
% sig to xy coordination system convert
% lonsig: Node number array, dimension [4,nodes]
% depsig: Layer depth array, demension [4,layers*nodes]

% Number of sigma layers.
fz = size(Mobj.siglayz, 2);
% Define the name of variables
% vartoplot = {'temperature','salinity'};
% Num = 24;

% Construct the data of plotting
[lonsig,depsig] = siglayz2xy(Mobj,crs_node);
% function [lonsig,depsig] = siglayz2xy(Mobj,crs_node)
% lonsig = [];
% depsig = [];
% lonsig(1:2,:) = repmat([1:size(crs_node,1)],2,fz);
% lonsig(3:4,:) = repmat([2:(size(crs_node,1)+1)],2,fz);
% depsig(1,:) = reshape(Mobj.siglayz(crs_node,:),[1,size(crs_node,1)*fz]);
% depsig(4,:) = reshape(Mobj.siglayz(crs_node,:),[1,size(crs_node,1)*fz]);
% depsig(2,(size(crs_node,1)+1):(size(crs_node,1)*fz)) = reshape(Mobj.siglayz(crs_node,1:(fz-1)),[1,size(crs_node,1)*(fz-1)]);
% depsig(3,(size(crs_node,1)+1):(size(crs_node,1)*fz)) = reshape(Mobj.siglayz(crs_node,1:(fz-1)),[1,size(crs_node,1)*(fz-1)]);

% Plot the animation
for i =1:length(vartoplot)
    fig = figure(Num+i);
    %clf
    fprintf([vartoplot{i},': cross section profile.','\n']);
    try
        minval = min(min(min(ncfile.(vartoplot{i})(crs_node,:,:))));
        maxval = max(max(max(ncfile.(vartoplot{i})(crs_node,:,:))));
    catch
        minval = min(min(ncfile.(vartoplot{i})(crs_node,:)));
        maxval = max(max(ncfile.(vartoplot{i})(crs_node,:)));
    end
    for tt =1:(timestep-1):length(ncfile.Times(:,1))
        try
            sudata = reshape(ncfile.(vartoplot{i})(crs_node,:,tt),[size(crs_node,1)*fz,1]);
        catch
            sudata = reshape(ncfile.(vartoplot{i})(crs_node,:,1),[size(crs_node,1)*fz,1]);
        end
        try
            DateString = ncfile.Times(tt,:);
        catch
            DateString = ncfile.Times(1,:);
        end
        plot_cross(Num+i, depsig, lonsig, sudata, (vartoplot{i}),'cross section',...
            minval,...
            maxval,...
            DateString);
        drawnow
        frame = getframe(fig);
        im{tt} = frame2im(frame);
        fprintf([vartoplot{i},': ',num2str(tt),' of time ',num2str(length(ncfile.Times)),'\n']);
    end
    % Save anamation files
    filename = ['plot_NC_cross_',vartoplot{i},'_',num2str(length(crs_node)),'.gif']; % Specify the output file name
    for tt =1:(timestep-1):length(ncfile.Times(:,1))
        [A,map] = rgb2ind(im{tt},256);
        if tt == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1/24);
        end
    end
end