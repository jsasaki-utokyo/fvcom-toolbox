function plot_nc_obc_TS(Num, Mobj, ncfile, timestep, vartoplot)
% Number of sigma layers.
fz = size(Mobj.siglayz, 2);
oNodes = [Mobj.read_obc_nodes{:}];
% Define the name of variables
% vartoplot = {'temperature','salinity'};
% Num = 24;
[lonsig,depsig] = siglayz2xy(Mobj,oNodes);

for i =1:length(vartoplot)
    fig = figure(Num+i);
    %clf
    fprintf([vartoplot{i},': boundary values.','\n']);
    minval = min(min(min(ncfile.(vartoplot{i})(oNodes,:,:))));
    maxval = max(max(max(ncfile.(vartoplot{i})(oNodes,:,:))));
    for tt = 1:(timestep-1):length(ncfile.Times)
        sudata = reshape(ncfile.(vartoplot{i})(oNodes,:,tt),[size(oNodes,1)*fz,1]);
        DateString = ncfile.Times(tt,:);
        plot_TS_cross(Num+i, depsig, lonsig, sudata, vartoplot{i},'open boundary',...
            minval,...
            maxval,...
            DateString);
        drawnow
        frame = getframe(fig);
        im{tt} = frame2im(frame);
        fprintf([vartoplot{i},': ',num2str(tt),' of time ',num2str(length(ncfile.Times)),'\n']);
    end
    % close;
    filename = ['plot_NC_obc_',vartoplot{i},'.gif']; % Specify the output file name
    for tt = 1:(timestep-1):length(ncfile.Times)
        [A,map] = rgb2ind(im{tt},256);
        if tt == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1/24/6);
        end
    end
end