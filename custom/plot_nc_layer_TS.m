function plot_nc_layer_TS(Num, ncfile, timestep, layers, vartoplot)
% sig to xy coordination system convert
% lonsig: Node number array, dimension [4,nodes]
% depsig: Layer depth array, demension [4,layers*nodes]

% Define the name of variables
% vartoplot = {'temperature','salinity'};
% Num = 24;

% For example
% vartoplot = {'temperature','salinity'};
% Num = 10;
% plot_nc_layer_TS(Num, ncfile, 4, [1], vartoplot);

% Plot the animation
for k = 1:length(vartoplot)
    for i = 1:length(layers)
        layer = layers(i);
        fig = figure(Num+k);
        %clf
        fprintf([vartoplot{k},': ',num2str(layer),' of layer ',num2str(length(ncfile.siglay(1,:))),'\n']);
        minval = min(min(min(ncfile.(vartoplot{k})(:,layer,:))));
        maxval = max(max(max(ncfile.(vartoplot{k})(:,layer,:))));
        for j = 1:(timestep-1):length(ncfile.Times(:,1))
            try
            plot_VALUE(Num+k, [ncfile.lon, ncfile.lat], ncfile.nv, ncfile.(vartoplot{k})(:,layer,j),...
                vartoplot{k}, ['layer No.',num2str(layer,'%02.2i')],...
                minval,...
                maxval,...
                [ncfile.Times(j,1:10),' ',ncfile.Times(j,12:19)]);
            catch
            plot_VALUE(Num+k, [ncfile.lon, ncfile.lat], ncfile.nv, ncfile.(vartoplot{k})(:,layer,j),...
                vartoplot{k}, ['layer No.',num2str(layer,'%02.2i')],...
                minval,...
                maxval,...
                [ncfile.Times(j,1:7)]);
            end
            drawnow
            frame = getframe(fig);
            im{j} = frame2im(frame);
            fprintf([vartoplot{k},': ',num2str(j),' of time ',num2str(length(ncfile.Times)),'\n']);
        end
        % close;
        filename = ['plot_NC_layer_',num2str(layer,'%02.2i'),'_',vartoplot{k},'.gif']; % Specify the output file name
        for j = 1:(timestep-1):length(ncfile.Times(:,1))
            [A,map] = rgb2ind(im{j},256);
            if j == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1/24);
            end
        end
    end
end