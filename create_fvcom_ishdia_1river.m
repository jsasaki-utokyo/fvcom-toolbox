%How to use?
%First;install TMD(TMD_Matlab_toolbox) from internet and install model data by visiting TPXO website
%Second; set path

%
% Which system am I using?
% CHANGE THESE TO YOUR OWN DIRECTOR
    basedir = 'C:/Users/Kota Ishida/';
      % Insert your mapped drive to \\store\login\your_name here
    addpath([basedir,'']);
    addpath([basedir,'Desktop/lab/fvcom-toolbox-main/fvcom_prepro/']);
    addpath([basedir,'Desktop/lab/fvcom-toolbox-main/utilities/']);
    addpath([basedir,'Desktop/lab/fvcom-toolbox-main/']);
    addpath([basedir,'Desktop/OceanMesh2D/work/']);
    addpath([basedir,'Desktop/lab/fvcom-toolbox-main/TMD_Matlab_Toolbox_v2.5/']);
    addpath([basedir,'Desktop/lab/fvcom-toolbox-main/TMD_Matlab_Toolbox_v2.5/TMD/']);
    addpath([basedir,'Desktop/lab/fvcom-toolbox-main/TMD_Matlab_Toolbox_v2.5/TMD/DATA/']);
    addpath([basedir,'Desktop/OceanMesh2D/datasets/RiverDischarge'])
    addpath([basedir,'Desktop/OceanMesh2D/work2/output/'])
    
% Base folder location - where do you want your forcing input files to end
% up?
inputConf.base = [basedir,'fvcominputs/test/'];
inputConf.outbase = [inputConf.base,'input/'];

%inputConf.AMM_folder = '\\store\projectsa\ISO_Modelling\POLCOMS\OPERATIONAL\OUTPUT\NETCDF\S12\AMM.hourly.';


% Which version of FVCOM are we using (for the forcing file formats)?
inputConf.FVCOM_version = '4.4.2';

% Case name for the model inputs and outputs
% Change this to whatever you want
inputConf.casename ='TokyoBay';
Mobj = read_sms_mesh('2dm','Tokyobay_v3f.2dm');
write_admesh_mesh(Mobj,'filename','out');
inputConf.grid = ['out.14'];
%Mobj = read_sms_mesh('2dm','Tokyo_bay_fix.2dm')
%write_admesh_mesh(Mobj,'filename','Tokyo_bay2')
%inputConf.grid = Mobj;
% output coordinates (FVCOM only likes cartesian at the moment)
inputConf.coordType = 'cartesian'; % 'spherical' or 'cartesian'
% input coordinates (what's my input bathy in?)
inputConf.coordInput = 'spherical'; % 'spherical' or 'cartesian'
% Input grid UTM Zone (if applicable)
inputConf.utmZone = {'54 N'};

% vertical coordinates type: sigma or hybrid
inputConf.verticalCoordType = 'sigma';

% Sponge layer parameters
inputConf.spongeRadius = -1; % in metres, or -1 for variable
inputConf.spongeCoeff = 0.001;

% z0 value in metres
inputConf.bedRoughness = 0.015; % or 0.015, 0.025 or 0.03 - Davies and Furnes (1980) shelf model

% Estimated velocity (m/s) and tidal range (m) for time step estimate
inputConf.estVel = 1.5;
inputConf.estRange = 2.0;

% Uniform temperature and salinity values
inputConf.temperature = 10;
inputConf.salinity = 35;

% Model time ([Y,M,D,h,m,s])
inputConf.modelYear = 2019;
inputConf.startDate = [inputConf.modelYear,8,10,18,00,00];
inputConf.endDate = [inputConf.modelYear,8,13,18,00,00];


% How many tidal constituents do we actually want to use at the model
% boundaries? Case sensitive (M2 != m2).
% Only relevant if using TPXO.
inputConf.tidalComponents = {'M2','S2','N2','K2','K1','O1','P1','Q1'};

% Give some names to the boundaries. This must match the number of node
% strings defined in SMS. Ideally, the order of the names should match the
% order in which the boundaries were made in SMS.
inputConf.boundaryNames = {'1'};

% Smooth the bathymetry if desired.
inputConf.smoothBathy = 'yes'; % 'yes' or 'no'.
if strcmpi(inputConf.smoothBathy, 'yes')
    % Set the smoothing factor and number of iterations (see smoothmesh).
    inputConf.smoothFactors = [0.5, 4]; % [factor, iterations]
end

%%%------------------------------------------------------------------------
%%%                      END OF INPUT CONFIGURATION
%%%------------------------------------------------------------------------

%% Prepare the data

% Convert times to Modified Julian Date
inputConf.startDateMJD = greg2mjulian(inputConf.startDate(1),inputConf.startDate(2),inputConf.startDate(3),inputConf.startDate(4),inputConf.startDate(5),inputConf.startDate(6));
inputConf.endDateMJD = greg2mjulian(inputConf.endDate(1),inputConf.endDate(2),inputConf.endDate(3),inputConf.endDate(4),inputConf.endDate(5),inputConf.endDate(6));
%inputConf.inputTimeTS = inputConf.startDateMJD:inputConf.dtTS:inputConf.endDateMJD;


% Read the input mesh and bathymetry. Also creates the data necessary for
% the Coriolis correction in FVCOM.
Mobj = read_grid_mesh('grid',inputConf.grid,...
    'coordinate',inputConf.coordType,'in_coord',inputConf.coordInput,...
    'project',true,'zone',inputConf.utmZone,'addCoriolis',true);

% Parse the open boundary nodes and add accordingly
% Add the sponge nodes
for i=1:size(Mobj.read_obc_nodes,2)
    nodeList = double(cell2mat(Mobj.read_obc_nodes(i)));
    Mobj = add_obc_nodes_list(Mobj,nodeList,inputConf.boundaryNames{i},1);
    if inputConf.spongeRadius < 0    % if we want a variable sponge radius
        if i==1
            % Create an array to store the radii
            Mobj.sponge_rad = zeros(size(Mobj.sponge_nodes));
        end
        % calculate the sponge radius
        spongeRadius = calc_sponge_radius(Mobj,nodeList);
        % Add the sponge nodes to the list
        Mobj = add_sponge_nodes_list(Mobj,nodeList,...
            [inputConf.boundaryNames{i},' sponge'],spongeRadius,...
            inputConf.spongeCoeff);
    else
        Mobj = add_sponge_nodes_list(Mobj,nodeList,...
            [inputConf.boundaryNames{i},' sponge'],inputConf.spongeRadius,...
            inputConf.spongeCoeff);
    end
    clear nodeList
end

clear i

%シグマ座標の取得
% Get the sigma depths in order to interpolate from the POLCOMS depths
if exist(fullfile(inputConf.outbase, 'sigma.dat'),'file')
    % If the sigma.dat file exists, read it
    Mobj = read_sigma(Mobj, fullfile(inputConf.base, 'input/sigma.dat'));
else
    % If we can't find the sigma.dat file, print an error message and
    % finish
    error(['sigma.dat not found. Please put your sigma.dat file into ',...
        fullfile(inputConf.outbase),' and try again.'])
end

% Do the bed roughness
Mobj.z0 = ones(1,Mobj.nElems)*inputConf.bedRoughness;


% Grid
write_FVCOM_grid(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_grd.dat']));

% Bathymetry
write_FVCOM_bath(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_dep.dat']));

% Coriolis
write_FVCOM_cor(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_cor.dat']));

% Open boundaries
write_FVCOM_obc(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_obc.dat']))

% Sponge file
write_FVCOM_sponge(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_spg.dat']))

% Bed roughness (constant or variable (see above))
write_FVCOM_z0(Mobj.z0,fullfile(inputConf.outbase,[inputConf.casename,'_z0.nc']),'bottom roughness');

%%


%%---------------------------------------------------
%                        tides
%---------------------------------------------------


%if strcmpi(inputConf.obcForcing, 'phase-amp')
    % Use elevation harmonics, not currents from TPXO
%    inputConf.extractType = 'z'; 
    % Need to cd to TPXO directory or it doesn't work
    % (Yes, this is inelegant but it's the easiest way for now)
%%    here = pwd; % store the current working directory to return later
 %   tpxo_dir = which('TMD');    % find the TPXO directory
 %   tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
 %   cd(tpxo_dir)    % go to TPXO directory
    % Location of the TMD model description file
    %inputConf.Model = [tpxo_dir,'DATA/Model_tpxo7.2'];
%elseif strcmpi(inputConf.obcForcing, 'z')
    % Need to cd to TPXO directory or it doesn't work
    % (Yes, this is inelegant but it's the easiest way for now)
 %%  tpxo_dir = which('TMD');    % find the TPXO directory
  %  tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
  %  cd(tpxo_dir)    % go to TPXO directory
     %Location of the TMD model description file
 %   inputConf.Model = ['Model_PO'];
%elseif strcmpi(inputConf.obcForcing, 'model-output')
    % Use NOC Operational Tide Surge Model output
%    inputConf.extractType = 'm';
%end
activate_Tide =1
if activate_Tide == 1
    inputConf.obcForcing = 'z'; 
    inputConf.datetide = 1/24
    inputConf.tidesMJD = inputConf.startDateMJD:inputConf.datetide:inputConf.endDateMJD;
    if strcmpi(inputConf.obcForcing, 'z')
        % Need to cd to TPXO directory or it doesn't work
        % (Yes, this is inelegant but it's the easiest way for now)
        here = pwd; % store the current working directory to return later
        tpxo_dir = which('TMD');    % find the TPXO directory
        tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
        cd(tpxo_dir)    % go to TPXO directory
        %Location of the TMD model description file
        inputConf.Model = ['Model_PO'];

        % Use tmd_tide_pred to predict surface elevations for a given time
        % range.
        
        % Add the tidal components to the Mobj.
        Mobj.Components = inputConf.tidalComponents;
        
        % Create a time series in MATLAB datenum format with ten minute
        % inputs
        inputConf.inputTimeTideZ = datenum(inputConf.startDate):1/144:datenum(inputConf.endDate);
        % Also do Modified Julian Day for the output to NetCDF
        inputConf.inputTimeTideZMJD = datenum(inputConf.startDateMJD):1/144:datenum(inputConf.endDateMJD);
        % Get the indices to use the tidal constituents defined in
        % inputConf.tidalComponents for TPXO (which requires a numerical
        % array of the constituents to be used). The order of the TPXO
        % constituents is M2, S2, N2, K2, K1, O1, P1, Q1, MF, MM, M4, MS4,
        % MN4.
        tpxoConsts = {'M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', ...
            'MF', 'MM', 'M4', 'MS4', 'MN4'};
        tIndUse = nan(length(Mobj.Components), 1);
        tInd = 1:length(tpxoConsts);
        for i=1:length(Mobj.Components)
            tPos = tInd(strcmp(Mobj.Components{i}, tpxoConsts));
            if ~isempty(tPos)
                tIndUse(i) = tPos;
            else
                warning('Supplied constituent (%s) is not present in the TPXO data', Mobj.Components{i}) %#ok<WNTAG>
            end
        end
        % Tidy up a bit
        clear c tpxoConsts tPos tInd
        % We can't just use tmd_tide_pred to do all the surface elevations
        % at once. Instead, the useful approaches are:
        %
        %   1. Time series at a single location
        %   2. Map of a given time step at all locations
        %
        % Since I'm likely to have many more time steps than locations,
        % it's probably best to do the time series at each location than
        % all the locations and a single time step.
        %
        % The order of the surface elevations in Mobj.surfaceElevation
        % should reflect the order of the open boundary node IDs as FVCOM
        % assumes they just map directly. So, rather than iterate through
        % each position, we need to get the position based on the list of
        % node IDs (Mobj.obc_nodes, without the zeros and in order
        % of each boundary).
        tmpObcNodes = Mobj.obc_nodes';
        ObcNodes = tmpObcNodes(tmpObcNodes~=0)';
        clear tmpObcNodes
    %         obc_lat = Mobj.lat(Mobj.obc_nodes(Mobj.obc_nodes~=0));
    %         obc_lon = Mobj.lon(Mobj.obc_nodes(Mobj.obc_nodes~=0));
        %Mobj.surfaceElevation = nan(size(ObcNodes,2), size(inputConf.inputTimeTideZ,2));
        %%for i=1:size(ObcNodes,2)
        %   currLon = Mobj.lon(ObcNodes(i));
        %   currLat = Mobj.lat(ObcNodes(i));
        %   fprintf('%,%',currrat,currLon);
        %end

    % for i=1:size(ObcNodes,2)
    %     currLon = Mobj.lon(ObcNodes(i));
    %     currLat = Mobj.lat(ObcNodes(i));
    %     fprintf('Position %i of %i (%.3f %.3f)... \n', i, size(ObcNodes,2), currLon, currLat)
    % end
    surfaceElevation = nan(Mobj.nObcNodes, size(inputConf.tidesMJD, 2), length(inputConf.boundaryNames));
    for i=1:length(inputConf.boundaryNames)
        for j=1:Mobj.nObcNodes
            % Get the current location (from the node ID)
            currLon = Mobj.lon(Mobj.obc_nodes(i,j));
            currLat = Mobj.lat(Mobj.obc_nodes(i,j));
            %if ftbverbose
                fprintf('Position %i of %i (%.3f %.3f)... \n', j, Mobj.nObcNodes, currLon, currLat);
            %end
            [surfaceElevation(j,:,i), ~] = tmd_tide_pred(inputConf.Model, ...
                inputConf.tidesMJD+678942.000000, currLat, currLon, 'z', tIndUse);
            if isnan(surfaceElevation(j,:))
                % Try the global model instead.
                [surfaceElevation(j,:,i), ~] = tmd_tide_pred(inputConf.Model, ...
                inputConf.tidesMJD, currLat, currLon, 'z', tIndUse);
            end
        end
    end
    for i=1:length(inputConf.boundaryNames)
        for j=1:Mobj.nObcNodes
            % Get the current location (from the node ID)
            currLon = Mobj.lon(Mobj.obc_nodes(i,j));
            currLat = Mobj.lat(Mobj.obc_nodes(i,j));
            %if ftbverbose
                fprintf('Position %i of %i (%.3f %.3f)... \n', j, Mobj.nObcNodes, currLon, currLat);
            %end
           % [surfaceElevation(j,:,i), ~] = [78,~];
        end
    end
    Mobj.surfaceElevation = surfaceElevation;
    % Tidy up some more
    clear i j tIndUse obc_lat obc_lon currLon currLat surfaceElevation
    cd(here);
    %save('varb/Mobj_01.mat','Mobj','-v7.3','-nocompression');

        %for i=1:size(ObcNodes,2)
            % Get the current location (from the node ID)
        %    currLon = Mobj.lon(ObcNodes(i));
        %   currLat = Mobj.lat(ObcNodes(i));
        %    fprintf('Position %i of %i (%.3f %.3f)... \n', i, size(ObcNodes,2), currLon, currLat)
        %    [Mobj.surfaceElevation(i,:), ~] = tmd_tide_pred(inputConf.Model, inputConf.inputTimeTideZ, currLat, currLon, 'z', tIndUse);
        %end
        % Tidy up a bit
        %clear tIndUse obc_lat obc_lon ObcNodes currLon currLat
        write_FVCOM_elevtide(Mobj, ...
        inputConf.tidesMJD,...
        fullfile(inputConf.outbase, [inputConf.casename, '_julian_obc.nc']),...
        'Model surface elevation boundary input',...
        'floattime', true,...
        'julian', true);
        fprintf('done.\n')
    else
        % No idea what has been supplied
        error('Unrecognised open boundary node forcing type. Choose phase')
    end
end





%%
%%%------------------------------------------------------------------------
%%%                     River discharge and output
%%%------------------------------------------------------------------------


Mobj.xc = nodes2elems(Mobj.x, Mobj);
Mobj.yc = nodes2elems(Mobj.y, Mobj);
Mobj.lonc = nodes2elems(Mobj.lon, Mobj);
Mobj.latc = nodes2elems(Mobj.lat, Mobj);



inputConf.riverForcing = 'FLUX';
% River information
% inputConf.river.infos = {...
%    'Tsurumigawa',...
%     'Sumidagawa',...
%     'Tamagawa',...
%     'Arakawa',...
%     'Edogawa',...
%     'Mamagawa',...
%     'Ebigawa',...
%     'Yorogawa',...
%     'Obitsugawa',...
%     'Koitogawa',...
%     'Muratagawa',...
%     'Hanamigawa'};
  inputConf.river.infos = {...
     %'Tsurumigawa',...
     %'Sumidagawa',...
     %'Tamagawa',...
     'Arakawa',...
     %'Edogawa',...
     %'Yorogawa',...
     %'Obitsugawa',...
     %'Koitogawa'
     };
% Location of river file
 inputConf.river.flux = [basedir,'Desktop/OceanMesh2D/datasets/RiverDischarge/riverflux_arakawa.csv'];
 inputConf.river.temp = [basedir,'Desktop/OceanMesh2D/datasets/RiverDischarge/rivertemp_arakawa.csv'];
 inputConf.river.salt = [basedir,'Desktop/OceanMesh2D/datasets/RiverDischarge/riversalt_arakawa.csv'];
 inputConf.river.location = [basedir,'Desktop/OceanMesh2D/datasets/RiverDischarge/riverloc_arakawa.csv'];

% Adjust river mouth location
% 139??55'57.84"	139??50'55.07"	139??46'29.73"	139??46'46.94"	139??40'53.63"	139??58'41.66"
%  35??41'56.29"	 35??38'36.31"	 35??38'49.41"	 35??31'45.11"	 35??28'25.39"	 35??40'52.25"
% New Edogawa river mouth location
% 139.872575 (139??52'21.27")
% 35.63695833(35??38'13.05")

% Give some names to the boundaries. This must match the number of node
% strings defined in SMS. Ideally, the order of the names should match the
% order in which the boundaries were made in SMS.
%inputConf.boundaryNames = {'pacific_ocean'};

% Stations
%inputConf.names = {...
%    'Tokyo',...
%    'Chiba',...
%    'Yokohamashinko',...
%    'Daini-Kaiho',...
%    'Yokosuka',...
%    'Kyurihamako',...
%     'Tidal-Station-1',...
%     'Tidal-Station-2',...
%     'Tidal-Station-3',...
%     'Tidal-Station-4',...
%    };
%inputConf.positions = [...
%    139.7700000000000,35.648888888888889;...
%    140.0455555555556,35.568055555555556;...
%    139.6441666666666,35.454166666666667;...
%    139.7433333333333,35.308611111111105;...
%    139.6513888888889,35.288055555555556;...
%    139.7208333333333,35.227777777777778;...
%     139.6912777777778,35.062205555555556;...
%     139.7207861111112,35.076750000000004;...
%     139.7504166666667,35.093661111111111;...
%     139.7735055555556,35.127727777777778;...
%    ];

% Make cross section courses
% 'Yoko_01': Tokyo-Chiba port;
% 'Yoko_02': Submarine high way
% 'Yoko_03': Inner bay boundary
% 'Yoko_04': Outer bay boundary
%inputConf.crosssection = {...
%   'Yoko_01',...
%    'Yoko_02',...
%    'Yoko_03',...
%    'Yoko_04',...
%    'Tate_01',...
%    'Tate_02',...
%    'Tate_03',...
%    };
%inputConf.endpoints = [...
%    139.8977,35.6166,140.0374,35.5166;...
%    139.8023,35.5112,139.9115,35.4378;...
%    139.7373,35.2656,139.7849,35.3121;...
%    139.6784,35.1397,139.7529,34.9784;...
%    140.0433,35.6303,139.7209,35.3020;...
%    139.7209,35.3020,139.7965,35.2395;...
%    139.7965,35.2395,139.7095,35.0545;...
%    ];
%inputConf.crosssection_ds = 0.01;

%[inputConf] = getCrossSectionPoints(inputConf);

% Generate coordinates
%UTMzone = regexpi(inputConf.utmZone,'\ ','split');
%for s = 1:length(inputConf.positions)
%    [inputConf.positions(s,3),inputConf.positions(s,4),~,~] = wgs2utm(...
%        inputConf.positions(s,2),inputConf.positions(s,1),...
%        str2double(char(UTMzone{1}(1))),char(UTMzone{1}(2)));
%end
%clear s UTMzone

% The accepted distance error of cal and real, in the cartesian coordinate 
% system, the dist unit is meter.
%inputConf.dist = 1100.00;


activate_River = 1
develop_mode = 1
if activate_River == 0
    disp('Skip creating river discharge files.')
else
if strcmpi(inputConf.riverForcing, 'FLUX')
    if develop_mode == 3
    	fprintf('Loading Model objet file...')
        load('varb/Mobj_04.mat');
        fprintf('Done!\n');
    else
        fprintf('Reading river forcing file...');
        Mobj = get_COSTUM_river_location(inputConf, Mobj);
        Mobj = get_COSTUM_river_variable(Mobj, ...
            {inputConf.river.flux,inputConf.river.temp,inputConf.river.salt}, ...
            {'flux','temp','salt'},...
            inputConf.river.infos,...
            'time', true);
        Mobj.nRivers = Mobj.river.number;
        inc = 1;
        Mobj.rivermouth = cell(1); % don't preallocate as we don't know how many we'll have
        for s = 1:Mobj.nRivers
            [node, ~] = find_nearest_pt(Mobj.river.location(s, 3), Mobj.river.location(s, 4), Mobj);
            [~, elem] = min(abs(sqrt((Mobj.xc - Mobj.river.location(s, 3)).^2 + Mobj.yc - Mobj.river.location(s, 4)).^2));
            Mobj.rivermouth{inc} = {inc, Mobj.river.location(s, 3), Mobj.river.location(s, 4), node, Mobj.h(node), Mobj.river.name(s), elem};
            riverList(s,1) = Mobj.rivermouth{inc}(1,4);
            Mobj.riverList = cell2mat(riverList);
            inc = inc + 1;
        end
        clear node elem inc s
        Mobj = add_river_nodes_list(Mobj,Mobj.riverList,Mobj.river.name,1);
        %save('varb/Mobj_04.mat','Mobj','-v7.3','-nocompression');
        fprintf('Done!\n');
        clear riverList;
    end
    % river file
    fprintf('Writing river forcing file...\n');
    write_FVCOM_river(fullfile(inputConf.outbase,...
        [inputConf.casename,'_river.nc']),...
         inputConf.river.infos,...
         Mobj.river.timeMJD,...
         Mobj.river.flux,...
         Mobj.river.temp + 2.0,...
         Mobj.river.salt,...
        'Tokyo Bay rivers',...
        'Model river boundary input');
    write_FVCOM_river_nml(Mobj, ...
        fullfile(inputConf.outbase,'RIVERS_NAMELIST.nml'), ...
        [inputConf.casename,'_river.nc'],...
        '''uniform''');
    fprintf('Done!\n');
end
end
%%
%clear ans tpxo_dir
fprintf('All done!\n')

