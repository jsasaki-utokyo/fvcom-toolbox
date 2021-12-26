%%
%%%------------------------------------------------------------------------
%%% Create FVCOM input files from a given OceanMesh2D "14" file
%%%     and Casename_sigma.dat file.
%%% Tested only for tidal forcing using tpxo9 and river forcing on Windows
%%%     for the purpose of storm surge simulation in Tokyo Bay.
%%% Place this script in fvcom-toolbox/Casenames/"Casename"/ along with
%%%     grid and sigma files, in which all create files will be written.
%%
%%
%%%------------------------------------------------------------------------
%%%                       DIRECTORY CONFIGURATION
%%%
%%%------------------------------------------------------------------------

clearvars; clc;

% If you are developing the code
% "1" means create file, download ncfile, and it takes time
%     This option is usually used when a new case is applied. 
%     All files will be generated in scratch.
% "2" means load downloaded files
%     This option is usually used when a new grid applied but in same
%     place, so that meteology and open boundary condition will not change.
%     For example: HYCOM_*.mat and FORCING_GWO.mat or FORCING_NCEP.mat
% "3" means load model object mat files
%     This option is usually used when running a same case, so that we dont
%     need to intorpolate hycom and forcing.
develop_mode = 2;

% If you wish to activate creating tidal forcing files:
% "1" means creating files.
% "0" means skipping creating files.
activate_tidal_forcing = 1;

% If you wish to activate creating HYCOM S & T forcing files:
% "1" means creating files.
% "0" means skipping creating files.
activate_HYCOM_S_T = 0;

% If you wish to activate creating River discharge files:
% "1" means creating files.
% "0" means skipping creating files.
activate_River = 1;

% Which system am I using?
if ismac    % On Mac
    basedir = '/Users/jsasaki/GitHub/';
    basedir2 = '/Volumes/jsasaki/';
    addpath([basedir,'fvcomtoolbox/']);
    addpath([basedir,'fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'fvcomtoolbox/utilities/']);
    addpath([basedir,'fvcomtoolbox/custom/']);
    addpath([basedir,'TMD/']);
    addpath([basedir,'TMD/FUNCTIONS/']);
    addpath([basedir,'TMD/DATA/']);
    %addpath([basedir2,'GitHub/paper_case/gwo/']);
    %addpath([basedir2,'GitHub/paper_case/RiverDischarge/']);
elseif isunix       % Unix?
    basedir = '/home/usr1/m70161a/';
    addpath([basedir,'Github/fvcomtoolbox/']);
    addpath([basedir,'Github/fvcomtoolbox/fvcom_prepro/']);
    addpath([basedir,'Github/fvcomtoolbox/utilities/']);
    addpath([basedir,'Github/fvcomtoolbox/custom/']);
    addpath([basedir,'Github/TMD/']);
    addpath([basedir,'Github/TMD/FUNCTIONS/']);
    addpath([basedir,'Github/TMD/DATA/']);
    addpath([basedir,'casesets/']);
    addpath([basedir,'casesets/gwo/']);
    addpath([basedir,'casesets/RiverDischarge/']);
elseif ispc     % Or Windows?
    basedir = 'D:\Github\';      
    addpath([basedir,'fvcom-toolbox\']);
    addpath([basedir,'fvcom-toolbox\fvcom_prepro\']);
    addpath([basedir,'fvcom-toolbox\utilities\']);
    addpath([basedir,'fvcom-toolbox\custom\']);
    addpath([basedir,'TMD_Matlab_Toolbox_v2.5\TMD\']);
    addpath([basedir,'TMD_Matlab_Toolbox_v2.5\TMD\FUNCTIONS\']);
    addpath([basedir,'TMD_Matlab_Toolbox_v2.5\TMD\DATA\']);
    addpath([basedir,'OceanMesh2D\datasets\RiverDischarge\']);
end

% Output directory: dat and nc files.
inputConf.outbase = pwd;

% Write out all the required files.
% Make the output directory if not exist.
if exist(inputConf.outbase, 'dir')~=7
    mkdir(inputConf.outbase)
end

% Working directory: grads folder and current folder
% inputConf.base = [basedir,'casesets/'];

% Write out diary file.
% Clear the diary file if it does exist.
if exist(fullfile(inputConf.outbase, 'diary'),'file')==2
    try
        diary off;
        delete('diary');
    catch
        delete('diary');
    end
end
diary on;

% Which version of FVCOM are we using (for the forcing file formats)?
%inputConf.FVCOM_version = '4.1';
%inputConf.FVCOM_version = '4.3.1';
inputConf.FVCOM_version = '4.4.2';
% Specify an OceanMesh2D "14" grid file
inputConf.grid = ['tb_futtsu.14'];

%%%------------------------------------------------------------------------
%%%                           Spatial stuff
%%%                           ????????????
%%%------------------------------------------------------------------------

% Case name for the model inputs and outputs
% Change this to whatever you want
inputConf.casename = 'tb_futtsu';

% output coordinates (FVCOM only likes cartesian at the moment)
inputConf.coordType = 'cartesian'; % 'spherical' or 'cartesian' 
% input coordinates (what's my input bathy in?)
inputConf.coordInput = 'spherical'; % 'spherical' or 'cartesian'

% Input grid UTM Zone (if applicable)
% See: https://upload.wikimedia.org/wikipedia/commons/e/ed/Utm-zones.jpg
% As Utm-zones indicated, here should be 54-S, but here the utmZone should
% be tmzone (UTM longitudinal zone) and utmhemi (UTM hemisphere as array of
% 'N' or 'S' characters)
inputConf.utmZone = {'54 N'};

% Option to smooth the bathymetry data.
inputConf.smoothBathy = 'yes'; % 'yes' or 'no'.
if strcmpi(inputConf.smoothBathy, 'yes')
    % Set the smoothing factor and number of iterations (see smoothmesh).
    inputConf.smoothFactors = [0.5, 4]; % [factor, iterations]
end

% vertical coordinates type: sigma or hybrid
inputConf.verticalCoordType = 'sigma';


%%%------------------------------------------------------------------------
%%%                     Time and Model constants
%%%                     ?????????????????????????????????
%%%------------------------------------------------------------------------

% Model time ([Y,M,D,h,m,s])
inputConf.modelYear = 2019;
inputConf.startDate = [inputConf.modelYear,10,04,18,00,00];
inputConf.endDate = [inputConf.modelYear,10,14,18,00,00];

% Convert times to Modified Julian Days
inputConf.startDateMJD = greg2mjulian(inputConf.startDate(1),...
    inputConf.startDate(2),inputConf.startDate(3),inputConf.startDate(4),...
    inputConf.startDate(5),inputConf.startDate(6));
inputConf.endDateMJD = greg2mjulian(inputConf.endDate(1),...
    inputConf.endDate(2),inputConf.endDate(3),inputConf.endDate(4),...
    inputConf.endDate(5),inputConf.endDate(6));

% The number of months in the period of data.
inputConf.mm = inputConf.startDate(2):inputConf.endDate(2);
if inputConf.mm == 1
    inputConf.dOffsets = [0, 4];
elseif inputConf.mm == 12
    inputConf.dOffsets = [2, 0];
else
    inputConf.dOffsets = [2, 4];
end

% Sponge layer parameters
inputConf.spongeRadius = -1; % in metres, or -1 for variable
inputConf.spongeCoeff = 0.001;

% z0 value in metres
inputConf.bedRoughness = 0.015; % or 0.015, 0.025 or 0.03 - Davies and Furnes (1980) shelf model

% Estimated velocity (m/s) and tidal range (m) for time step estimate
inputConf.estVel = 1.5;
inputConf.estRange = 2.0;


%%%------------------------------------------------------------------------
%%%                       Forcing and stuff
%%%                       ?????????????????????
%%%------------------------------------------------------------------------

% Model time type ('non-julian' and 'julian')
inputConf.datetype = 'julian';
% Increment used for tide (days)
inputConf.datetide = 1/24;
% Model time is decided dependent on the datetype.
if strcmpi(inputConf.datetype,'non-julian') 
    inputConf.modelTime = [...
        inputConf.startDateMJD, ...
        inputConf.endDateMJD];
elseif strcmpi(inputConf.datetype,'julian')
    inputConf.modelTime = [...
        inputConf.startDateMJD - inputConf.dOffsets(1), ...
        inputConf.endDateMJD + inputConf.dOffsets(2)];
end
% Open boundary forcing nodal forcing type (drived by OTPS).
inputConf.obcForcing = 'z'; 
inputConf.tidesMJD = inputConf.startDateMJD:inputConf.datetide:inputConf.endDateMJD;

% Increment used for open boundary ST (days)
inputConf.dateobs = 1/24;
% Open boundary temperatures and salinities (string for source or number for constant).
% The data is avaliable from 1992-10-02 00:00:00
inputConf.obc_temp = 'HYCOM';
inputConf.obc_salt = 'HYCOM';
inputConf.obc_u = 'NONE';
inputConf.obc_v = 'NONE';
inputConf.obctsMJD = [inputConf.startDateMJD, inputConf.endDateMJD + 1];
if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
  inputConf.Nested_type = 1;% [1, 2 == direct nesting, 3 == weighted]
  inputConf.levels = 1;% number of boundary bands to use for weighted option
  inputConf.power = 0;% [0 is linear, anything else is 1/conf.levels.^conf.power]
end

% Increment used for surface forcing (days)
inputConf.dateforcing = 1/24;
% The surface forcing from NCEP
% The data is avaliable from 1948 - present
inputConf.doForcing = 'GWO';
if strcmpi(inputConf.doForcing, 'GWO')
    inputConf.forceMJD = inputConf.startDateMJD:inputConf.dateforcing:inputConf.endDateMJD;
elseif strcmpi(inputConf.doForcing, 'NCEP')
    inputConf.forceMJD = [inputConf.startDateMJD, inputConf.endDateMJD];
end

% The river forcing from flux, salinity and temperature
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
     'Tsurumigawa',...
     'Sumidagawa',...
     'Tamagawa',...
     'Arakawa',...
     'Edogawa',...
     'Yorogawa',...
     'Obitsugawa',...
     'Koitogawa'};
% Location of river file
 inputConf.river.flux = [basedir,'OceanMesh2D/datasets/RiverDischarge/river_flux.csv'];
 inputConf.river.temp = [basedir,'OceanMesh2D/datasets/RiverDischarge/river_temp.csv'];
 inputConf.river.salt = [basedir,'OceanMesh2D/datasets/RiverDischarge/river_salt.csv'];
 inputConf.river.location = [basedir,'OceanMesh2D/datasets/RiverDischarge/river_location.csv'];

% Adjust river mouth location
% 139??55'57.84"	139??50'55.07"	139??46'29.73"	139??46'46.94"	139??40'53.63"	139??58'41.66"
%  35??41'56.29"	 35??38'36.31"	 35??38'49.41"	 35??31'45.11"	 35??28'25.39"	 35??40'52.25"
% New Edogawa river mouth location
% 139.872575 (139??52'21.27")
% 35.63695833(35??38'13.05")

% Give some names to the boundaries. This must match the number of node
% strings defined in SMS. Ideally, the order of the names should match the
% order in which the boundaries were made in SMS.
inputConf.boundaryNames = {'pacific_ocean'};

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

%%
%%%------------------------------------------------------------------------
%%%                      Model mesh generation
%%%                      ????????????
%%%------------------------------------------------------------------------

% Read the input mesh and bathymetry. Also creates the data necessary for
% the Coriolis correction in FVCOM.
Mobj = read_grid_mesh(...
    'grid',inputConf.grid,...
    'coordinate',inputConf.coordType,...
    'in_coord',inputConf.coordInput,...
    'project',true,...
    'zone',inputConf.utmZone,...
    'addCoriolis',true);

% Check the consistency of open boundaries
disp('Check the consistency')
if length(Mobj.read_obc_nodes) ~= length(inputConf.boundaryNames)
    length(Mobj.read_obc_nodes)
    disp('ERROR: Mismatch in the number of open boundaries in inputConf.boundaryNames')
    return
end

% Smooth the bathymetry if desired.
if strcmpi(inputConf.smoothBathy, 'yes')
    Mobj = setup_metrics(Mobj);
    % Backup Mobj.h as Mobj.h_backup
    Mobj.h_backup_1 = Mobj.h;
    Mobj.h = smoothfield(Mobj.h, Mobj, ...
        inputConf.smoothFactors(1), inputConf.smoothFactors(2));
    Mobj.h_backup_2 = Mobj.h;
    Mobj.h = smoothfield(Mobj.h, Mobj, ...
        inputConf.smoothFactors(1), inputConf.smoothFactors(2));
    % smoothfield2 is really inappropriate for bathymetry data.
    % Mobj.h = smoothfield2(Mobj.h,Mobj,inputconf.smoothFactors(2));
    Mobj.h_delta = Mobj.h - Mobj.h_backup_1;
    plotMesh(01, [Mobj.lon, Mobj.lat], Mobj.tri, Mobj.h_delta);
end

% % Parse the open boundary nodes and add accordingly
% % Add the sponge nodes
% x_adj = Mobj.x;
% y_adj = Mobj.y;
% for i=1:size(Mobj.read_obc_nodes,2)
%     % For each of the open boundaries, find the points inside the
%     % boundary and move the closest existing point to that location.
%     % This should make the boundaries more orthogonal (as suggested by
%     % the FVCOM manual in Figure 6.2 on page 77 [as of 2013-03-11]).
%     % x= x3;
%     % y= y3;
%     % node_ids=Mobj.read_obc_nodes{1};
%     [x_adj, y_adj] = fix_inside_boundary(x_adj, y_adj, Mobj.read_obc_nodes{i});    
% end
% Mobj.x = x_adj;
% Mobj.y = y_adj;

% Parse the open boundary nodes and add accordingly
% Add the sponge nodes
for i=1:size(Mobj.read_obc_nodes,2)
    nodeList = double(cell2mat(Mobj.read_obc_nodes(i)));
    Mobj = add_obc_nodes_list(Mobj,nodeList,inputConf.boundaryNames{i},1,1);
    % find inside open bounday
    % [~,~,temp_list] = fix_inside_boundary(Mobj.x, Mobj.y, Mobj.read_obc_nodes{i});
    % disp(temp_list)
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
    clear nodeList spongeRadius
end
clear i
%save('varb/temp_list.mat','temp_list','-v7.3','-nocompression');
% plot points inside open boundary
%{
clear i
x = Mobj.x;
y = Mobj.y;
plot(x(temp_list), y(temp_list), 'go');
%}

%
if strcmpi(inputConf.verticalCoordType,'sigma')
    % Get the sigma depths in order to interpolate from the depths
    if exist(fullfile(inputConf.outbase, [inputConf.casename,'_sigma.dat']),'file')
        % If the sigma.dat file exists, read it
        Mobj = read_sigma(Mobj, fullfile(inputConf.outbase, [inputConf.casename,'_sigma.dat']));
    else
        % If we can't find the sigma.dat file, print an error message and finish
        error(['Casename_sigma.dat not found. Please put your sigma.dat file into ',...
            fullfile(inputConf.outbase),' and try again.'])
    end
elseif strcmpi(inputConf.verticalCoordType,'hybrid')
    inputConf.hybrid.sigma_file = [inputConf.casename,'_hybrid.dat'];
    inputConf.hybrid.nlev = 11;                % number of vertical levels (layers + 1)
    inputConf.hybrid.H0 = 100;                 % transition depth of the hybrid coordinates
    inputConf.hybrid.KU = 2;                   % number of layers in the DU water column
    inputConf.hybrid.KL = 1;                   % number of layers in the DL water column
    inputConf.hybrid.DU = 20;                  % upper water boundary thickness (metres)
    inputConf.hybrid.DL = 10;                  % lower water boundary thickness (metres)
    Mobj = hybrid_coordinate(inputConf.hybrid, Mobj);
end

% Do the bed roughness
Mobj.z0 = ones(1,Mobj.nElems)*inputConf.bedRoughness;

% Generate center point coordination of elements
% Do the stations list
Mobj.xc = nodes2elems(Mobj.x, Mobj);
Mobj.yc = nodes2elems(Mobj.y, Mobj);
Mobj.lonc = nodes2elems(Mobj.lon, Mobj);
Mobj.latc = nodes2elems(Mobj.lat, Mobj);

% Add station
% Positions=inputConf.positions;Names=inputConf.names;Dist=inputConf.dist;plotFig=1;
%Mobj = add_stations_list(Mobj,inputConf.positions,inputConf.names,inputConf.dist,1);

% Estimate model time step. Supply estimated velocity (m/s) and tidal range
% (m) after the mesh object.
Mobj = estimate_ts(Mobj,inputConf.estVel,inputConf.estRange);
fprintf('Estimated time step:\t%.2f\n',min(Mobj.ts));

%%
%%%------------------------------------------------------------------------
%%%                     Output of basic configurations
%%%                     ????????????
%%%------------------------------------------------------------------------

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
% Time series wave stations
%write_FVCOM_stations(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_station.dat']));
% 2dm file for gis
%write_SMS_2dm(fullfile(inputConf.outbase,[inputConf.casename,'_grd.2dm']), Mobj.tri, Mobj.x, Mobj.y, Mobj.h);

% Save Model object file
if exist('varb', 'dir')~=7
    mkdir('varb')
end
save('varb/Mobj_00.mat','Mobj','-v7.3','-nocompression');

 
%%
%%%------------------------------------------------------------------------
%%%                    Additional forcing: Tides
%%%                    ??????????????? 
%%%------------------------------------------------------------------------
if activate_tidal_forcing == 0
    disp('Skip creating tidal forcing files (activate_tidal_forcing = 0)')
else
tic
% Open boundary nodal forcing type
%   'z' for predicted surface elevation
%   'phase-amp' for amplitudes and phases
%   'model-output' for Tidal Model Driver output
fprintf('Calculating open boundary forcing data from OTPS...\n')
% Need to cd to TPXO directory or it doesn't work
inputConf.extractType = 'z'; 
% (Yes, this is inelegant but it's the easiest way for now)
here = pwd; % store the current working directory to return later
tpxo_dir = which('TMD');    % find the TPXO directory
tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
disp(tpxo_dir)
if strcmpi(inputConf.obcForcing, 'z')
    cd(tpxo_dir)    % go to TPXO directory
    % Location of the TMD model description file
    %inputConf.Model = [tpxo_dir,'DATA/Model_OhS'];
    inputConf.Model = [tpxo_dir,'DATA/Model_tb_futtsu'];
elseif strcmpi(inputConf.obcForcing, 'phase-amp')
    cd(tpxo_dir)    % go to TPXO directory
    % Location of the TMD model description file
    %inputConf.Model = [tpxo_dir,'DATA/Model_OhS'];
    inputConf.Model = [tpxo_dir,'DATA/Model_tb_futtsu'];
elseif strcmpi(inputConf.obcForcing, 'otps')
    fid=fopen([tpxo_dir,'LAT_LON/lat_lon_1st'],'w');
    fprintf('Rewriting %s...\n','lat_lon_1st');
    fprintf(fid,"   lat       lon             yy   mm   dd   hh   mi  sec  dt(min) TSLength\n");
    for i = 1:Mobj.nObcNodes
        fprintf(fid," %10.4f %10.4f %8d %4d %4d %4d %4d %4d %4d %10d\n",...
            Mobj.lat(i), Mobj.lon(i),...
            inputConf.startDate(1), inputConf.startDate(2),...
            inputConf.startDate(3), inputConf.startDate(4),...
            inputConf.startDate(5), inputConf.startDate(6),...
            inputConf.datetide*24*60, 1/inputConf.datetide*(inputConf.endDateMJD-inputConf.startDateMJD)+1);
    end
    fclose(fid);
    fprintf('Finishing rewriting %s...\n','lat_lon_1st');
end

% How many tidal constituents do we actually want to use at the model
% boundaries? Case sensitive (M2 != m2).
% Only relevant if using TPXO.
% inputConf.tidalComponents = {'M2','S2','N2','K2','K1','O1','P1','Q1'};
inputConf.tidalComponents = {'M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1','MF','MM','M4','MS4','MN4','2N2','S1'};
% clear tpxo_dir

%%%------------------------------------------------------------------------
%%%                  Tides and its (non)julian output
%%%------------------------------------------------------------------------

% Generate a surface elevation time series for open boundary forcing.
if strcmpi(inputConf.obcForcing, 'z')
    if develop_mode == 3
        cd(here);
    	fprintf('Loading Model objet file...\n');
        load('varb/Mobj_01.mat');
        fprintf('Done!\n');
    else
        % Use tmd_tide_pred to predict surface elevations for a given time range.
        % Add the tidal components to the Mobj.
        % Mobj.Components = conf.obc_tides.components;
        Mobj.Components = inputConf.tidalComponents;
        % Get the indices to use the tidal constituents defined in
        % conf.obc_tides.components for TPXO (which requires a
        % numerical array of the constituents to be used). The order of the
        % TPXO constituents is M2, S2, N2, K2, K1, O1, P1, Q1, MF, MM, M4,
        % MS4, MN4.
        % tpxoConsts = {'M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1'};
        tpxoConsts = {'M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1','MF','MM','M4','MS4','MN4','2N2','S1'};

        tIndUse = nan(length(Mobj.Components), 1);
        tInd = 1:length(tpxoConsts);
        for i=1:length(Mobj.Components)
            tPos = tInd(strcmp(Mobj.Components{i}, tpxoConsts));
            if ~isempty(tPos)
                tIndUse(i) = tPos;
            else
                warning('Supplied constituent (%s) is not present in the TPXO data', Mobj.Components{i}); %#ok<WNTAG>
            end
        end
        % Tidy up a bit
        clear i c tpxoConsts tPos tInd
        tIndUse = tIndUse(~isnan(tIndUse));
        % We can't just use tmd_tide_pred to do all the surface elevations
        % at once. Instead, the useful approaches are:
        %   1. Time series at a single location
        %   2. Map of a given time step at all locations
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
        Mobj.surfaceElevation = surfaceElevation;
        % Tidy up some more
        clear i j tIndUse obc_lat obc_lon currLon currLat surfaceElevation
        cd(here);
        save('varb/Mobj_01.mat','Mobj','-v7.3','-nocompression');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEST CODE
        % t.Model = '/home/usr0/n70110d/github/fvcomtoolbox/tmd/DATA/Model_OhS'
        % t.SDtime = [datenum([2019,05,14,15,00,00]):1/24/60*5:datenum([2019,05,14,15,00,00])+1]
        % t.lat = 35.5681
        % t.lon = 140.0456
        % t.ptype = 'z'
        % t.Cid = 
        % [t.TS,t.conList]=tmd_tide_pred(t.Model,t.SDtime,t.lat,t.lon,t.ptype)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
elseif strcmpi(inputConf.obcForcing, 'otps')
    if develop_mode == 3
    	fprintf('Loading Model objet file...\n');
        load('varb/Mobj_01.mat');
        fprintf('Done!\n');
    else
        surfaceElevation = nan(Mobj.nObcNodes, size(inputConf.tidesMJD, 2), length(inputConf.boundaryNames));
        for i=1:length(inputConf.boundaryNames)
            for j=1:Mobj.nObcNodes
                % Get variables from the generated mat format files
                filename = [tpxo_dir,'OUT/1st_',mat2str(j),'.mat'];
                fprintf('Getting %s from TMD\n',['1st_',mat2str(j),'.mat']);
                variable = {'SerialDay','TimeSeries'};
                load(filename,variable{:});
                [~,n] = find(SerialDay - 678942 == inputConf.tidesMJD(1));
                for k=1:length(inputConf.tidesMJD)
                    surfaceElevation(j,k,i) = TimeSeries(1,n+k-1);
                end
            end
        end
        Mobj.surfaceElevation = surfaceElevation;
        % Tidy up some more
        clear i j k m n filename variable TimeSeries SerialDay surfaceElevation
        save('varb/Mobj_01.mat','Mobj','-v7.3','-nocompression');
    end
elseif strcmpi(inputConf.obcForcing,'phase-amp')
    if develop_mode == 3
        cd(here);
    	fprintf('Loading Model objet file...\n');
        load('varb/Mobj_01.mat');
        fprintf('Done!\n');
    else
        % Boundary conditions from TPXO (for spectral tides or predicted
        % surface elevations)
        % Put the input list into the mesh object.
        Mobj.Components = inputConf.tidalComponents;
        % Set up the tidal struct. This contains the relevant information for up to
        % eight constituents, ordered as period, beta love number and equilibrium
        % amplitude.
        %                    period    beta    eq. amp.
        %                      (s)      (?)       (m)
        tideComponents.M2 = [44714.16, 0.693, 0.242334];
        tideComponents.S2 = [43200.00, 0.693, 0.112841];
        tideComponents.N2 = [45570.24, 0.693, 0.046398];
        tideComponents.K2 = [43082.28, 0.693, 0.030704];
        tideComponents.K1 = [86163.84, 0.736, 0.141565];
        tideComponents.O1 = [92949.84, 0.695, 0.100514];
        tideComponents.P1 = [86637.24, 0.706, 0.046843];
        tideComponents.Q1 = [96726.24, 0.695, 0.019256];
        %tideComponents.Mf = [1180260,  ?????, 0.041742];
        %tideComponents.Mm = [2380716,  ?????, 0.022026];
        %tideComponents.Ssa = [15778980, ????, 0.019446];

        % Extract the values for each tidal component into Mobj.period_obc,
        % Mobj.beta_love and Mobj.equilibrium_amp.
        for c=1:size(Mobj.Components,2)
            Mobj.period_obc(c) = tideComponents.(Mobj.Components{c})(1);
            Mobj.beta_love(c) = tideComponents.(Mobj.Components{c})(2);
            Mobj.equilibrium_amp(c) = tideComponents.(Mobj.Components{c})(3);
        end
        clear c
        % Provide amplitude and phase data for the boundary nodes. Use the TMD
        % function tmd_extract_HC.m to get harmonic constants at the boundary
        % nodes.
        amp=cell(1,Mobj.nObs);
        Gph=cell(1,Mobj.nObs);
        Depth=cell(1,Mobj.nObs);
        constList=cell(1,Mobj.nObs);
        for i=1:length(inputConf.boundaryNames)
            % It is possible to specify the indices of the constituents of interest
            % when calling tmd_extract_HC, but it requires knowing the order
            % they're stored in the file. Easier for me to extract the constituents
            % of interest separately. This makes it a bit slower (having to
            % interpolate all the constituents is slower than a select few), but
            % it's a bit easier to code up.
            if Mobj.have_lonlat
            [amp{i},Gph{i},Depth{i},constList{i}] = tmd_extract_HC(inputConf.Model,Mobj.lat(Mobj.read_obc_nodes{i}),Mobj.lon(Mobj.read_obc_nodes{i}),inputConf.extractType);
            else
                % Need to convert XY to latlon.
                try % to use the handy file exchange utm2deg function
                    % Make cell array of all the zones because utm2deg is a bit
                    % inflexible in that regard (size of utmZones must equal size
                    % of x and y).
                    % This is somewhat redundant now that the lat/long is added
                    % when generating the Coriolis values, but it's still
                    % worthwhile keeping it here just in case. No harm should
                    % come of it being here anyway.
                    utmZones=cellfun(@(x) repmat(x,length(Mobj.x(Mobj.read_obc_nodes{i})),1),inputConf.utmZone,'uni',false);
                    [tmpLat,tmpLon] = utm2deg(Mobj.x(Mobj.read_obc_nodes{i}),Mobj.y(Mobj.read_obc_nodes{i}),utmZones{1});
                    % Get the tidal data
                    [amp{i},Gph{i},Depth{i},constList{i}] = tmd_extract_HC(inputConf.Model,tmpLat,tmpLon,inputConf.extractType);
                catch %#ok<CTCH>
                    error('Can''t convert X/Y positions to lat/long, so can''t extract data from the TPXO data. Consider adding utm2deg to your PATH.')
                end
            end
            for j=1:numel(Mobj.Components)
                fprintf('Extracting %s... ',Mobj.Components{j});
                posIdx = strmatch(lower(Mobj.Components{j}),constList{i}); %#ok<MATCH2>
                Mobj.amp_obc{i}(j,:) = amp{i}(posIdx,:);
                Mobj.phase_obc{i}(j,:) = Gph{i}(posIdx,:); % Greenwich phase
                fprintf('done!\n');
            end
        end
        clear posIdx amp Gph Depth constList
        % Find NaNs in the boundaries
        for i=1:Mobj.nObs
            brokenBoundary=i;
            nanIdx = Mobj.read_obc_nodes{brokenBoundary}(isnan(Mobj.phase_obc{brokenBoundary}(1,:)));
            nanLon = Mobj.lon(nanIdx);
            nanLat = Mobj.lat(nanIdx);
            inputConf.doFig=0;
            if max(nanLon)-min(nanLon)==0
                minPos = min(nanLat);
                maxPos = max(nanLat);
                inputConf.doFig=1;
            elseif max(nanLat)-min(nanLat)==0
                minPos = min(nanLon);
                maxPos = max(nanLon);
                inputConf.doFig=1;
            elseif isempty(nanIdx)
                fprintf('No NaNs in %s boundary.\n',inputConf.boundaryNames{i});
                clear nanLon nanLat nanIdx
            else
                error('Boundaries are not linear. Won''t plot %s boundary NaNs',inputConf.boundaryNames{i});
            end
            if inputConf.doFig
                figure
                patch('Vertices',[Mobj.lon,Mobj.lat],'Faces',Mobj.tri,...
                    'Cdata',Mobj.h,'edgecolor','k','facecolor','interp');
                hold on;
                plot(Mobj.lon(nanIdx),Mobj.lat(nanIdx),'wo','LineWidth',3,'MarkerSize',12);
                plot(Mobj.lon(nanIdx),Mobj.lat(nanIdx),'ko','LineWidth',3,'MarkerSize',8);
                axis('equal','tight');
            end
        end
        save('varb/Mobj_01.mat','Mobj','-v7.3','-nocompression');
    end
end
clear tpxo_dir;

if strcmpi(inputConf.obcForcing, 'z')
    % Write out the TPXO predicted surface elevation.
    write_FVCOM_elevtide(Mobj, ...
        inputConf.tidesMJD,...
        fullfile(inputConf.outbase, [inputConf.casename, '_julian_obc.nc']),...
        'Model surface elevation boundary input',...
        'floattime', true,...
        'julian', true);
elseif strcmpi(inputConf.obcForcing, 'otps')
    % Write out the TPXO predicted surface elevation.
    write_FVCOM_elevtide(Mobj, ...
        inputConf.tidesMJD,...
        fullfile(inputConf.outbase, [inputConf.casename, '_julian_obc.nc']),...
        'Model surface elevation boundary input',...
        'floattime', true,...
        'julian', true);
elseif strcmpi(inputConf.obcForcing,'phase-amp')  
    % Write out the TPXO predicted spectral tide.
    set_spectide(Mobj,...
        numel(Mobj.Components),...
        fullfile(inputConf.outbase,[inputConf.casename,'_non_julian_obc.nc']),...
        'TPXO spectral tidal boundary input');
end
fprintf('Tidal forcing working time: %.2f minutes\n', toc / 60);

end

if activate_HYCOM_S_T == 0
    disp('Skip HYCOM S & T forcing')
else
%%
%%%------------------------------------------------------------------------
%%%                    HYCOM S&T forcing and staff
%%%------------------------------------------------------------------------
tic
% Now we need some boundary temperature and salinity conditions.
if strcmpi('HYCOM', {inputConf.obc_temp, inputConf.obc_salt})
    % Use HYCOM data for the boundary forcing.
    % Offset the times to give us a bit of wiggle room.
    if develop_mode == 3
        fprintf('Loading Model objet file...\n')
        load('varb/Mobj_02.mat');
        fprintf('Done!\n');
    else
        if develop_mode == 1
            % modelTime = inputConf.obctsMJD;varlist={'temperature', 'salinity'}
            fprintf('Downloading daliy open boundary S&T forcing from HYCOM...\n');
            % T
            for i = 1:10
                try
                    hycom_t = get_HYCOM_forcing(Mobj, inputConf.obctsMJD, {'temperature'}); 
                    break;  % Break out of the i-loop on success
                catch ME
                    disp(ME);
                    fprintf('Retrying...\n');
                end
            end
            save(['../00_data/hycom_t','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat'],...
                'hycom_t','-v7.3','-nocompression');
            % S
            for i = 1:10
                try
                    hycom_s = get_HYCOM_forcing(Mobj, inputConf.obctsMJD, {'salinity'}); 
                    break;  % Break out of the i-loop on success
                catch ME
                    disp(ME);
                    fprintf('Retrying...\n');
                end
            end
            save(['../00_data/hycom_s','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat'],...
                'hycom_s','-v7.3','-nocompression');
            fprintf('Downloading daliy open boundary meanflow from HYCOM...\n');
            if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
                fprintf('Writing daliy open boundary meanflow file.\n')
                % u
                for i = 1:10
                    try
                        hycom_u = get_HYCOM_forcing(Mobj, inputConf.obctsMJD, {'u'}); 
                        break;  % Break out of the i-loop on success
                    catch ME
                        disp(ME);
                        fprintf('Retrying...\n');
                    end
                end
                save(['../00_data/hycom_u','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat'],...
                    'hycom_u','-v7.3','-nocompression');
                % v
                for i = 1:10
                    try
                        hycom_v = get_HYCOM_forcing(Mobj, inputConf.obctsMJD, {'v'}); 
                        break;  % Break out of the i-loop on success
                    catch ME
                        disp(ME);
                        fprintf('Retrying...\n');
                    end
                end
                save(['../00_data/hycom_v','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat'],...
                    'hycom_v','-v7.3','-nocompression');
            end
            fprintf('Downloading daliy open boundary S&T forcing from HYCOM...Done!\n');
        elseif develop_mode == 2
            fprintf('Loading daliy open boundary S&T forcing from local HYCOM database...\n')
            load(['../00_data/hycom_t','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat']);
            load(['../00_data/hycom_s','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat']);
            fprintf('Downloading daliy open boundary S&T forcing from HYCOM...Done!\n');
            if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
                load(['../00_data/hycom_u','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat']);
                load(['../00_data/hycom_v','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat']);
            end
        end
        % Interpolate the 4D HYCOM data on the FVCOM vertical grid at the open boundaries.
        Mobj = get_HYCOM_tsobc(Mobj, hycom_t, {'temperature'});
        Mobj = get_HYCOM_tsobc(Mobj, hycom_s, {'salinity'});
        if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
            % finding nesting region shows a nand of elements,
            % that is not we wanted, we need obc faces, instead.
            % Nested = find_nesting_region(inputConf, Mobj);
            % Mobj.read_obc_elems = Nested.read_obc_elems;
            Mobj = find_boundary_elements(Mobj);
            Mobj = get_HYCOM_tsobc(Mobj, hycom_u, {'u'});
            Mobj = get_HYCOM_tsobc(Mobj, hycom_v, {'v'});
        end
        clear hycom_*
        % backup daliy data
        Mobj.backup_temp = Mobj.temperature; 
        Mobj.backup_salt = Mobj.salt;
        Mobj.backup_tstm = Mobj.ts_times;
        if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
            Mobj.backup_mflu = Mobj.u;
            Mobj.backup_mflv = Mobj.v;
            Mobj.backup_mftm = Mobj.mf_times;
        end
        % recover daliy data
        %{
        Mobj.temperature = Mobj.backup_temp; Mobj.salt = Mobj.backup_salt;
        Mobj.u = Mobj.backup_mflu;           Mobj.v = Mobj.backup_mflv;
        Mobj.ts_times = Mobj.backup_tstm;    Mobj.mf_times = Mobj.backup_mftm;
        %}
        % Interpolate the 4D HYCOM data on the hourly time series
        Mobj = get_HYCOM_series(Mobj, inputConf.dateobs, 'temperature',true);
        Mobj = get_HYCOM_series(Mobj, inputConf.dateobs, 'salinity',true);
        if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
            Mobj = get_HYCOM_series(Mobj, inputConf.dateobs, 'u',true);
            Mobj = get_HYCOM_series(Mobj, inputConf.dateobs, 'v',true);
            % Generate boundary mean flow values at the centroids of the open boundary elements. 
            % For the various output files, we need both u and v as well as depth averaged velocity 
            % at each boundary node position. 
            % An example: get_POLCOMS_meanflow uses PML POLCOMS-ERSEM NetCDF files 
            % to interpolate mean flow to the FVCOM open boundary elements and vertical grid. 
            % The velocity data should be saved in Mobj.velocity of size [nElements, nTime] 
            % and the u and v components in Mobj.meanflow_u and Mobj.meanflow_v as arrays 
            % of size [nElements, nSiglay, nTime].
            % Mobj.meanflow = zeros(numel(Mobj.read_obc_elements{1}), numel(Mobj.mf_times));
            % Stick the values in the mesh structure.
            Mobj.meanflow_u = Mobj.u;
            Mobj.meanflow_v = Mobj.v;
            % Now we have the 3D arrays, create depth averaged velocities too
            Mobj.meanflow_ubar = squeeze(mean(Mobj.meanflow_u, 2));
            Mobj.meanflow_vbar = squeeze(mean(Mobj.meanflow_v, 2));
            % Depth averaged velocity
            Mobj.velocity = squeeze(mean(sqrt(Mobj.meanflow_u.^2 + Mobj.meanflow_v.^2), 2));
        end
        save('varb/Mobj_02.mat','Mobj','-v7.3','-nocompression');
    end
elseif strcmpi('FRA-JCOPE', {inputConf.obc_temp, inputConf.obc_salt})
    % [data,header]=read_grads(file_name,var_name,varargin)
    % 
    % file_name = ['C:\Users\Yulong WANG\Documents\GitHub\jcope-convert\fra_jcope\el.ctl'];
    % var_name = ['all'];
    % [data,header]=read_grads('C:\Users\Yulong WANG\Documents\GitHub\jcope-convert\fra_jcope\t.ctl','all'); 
end
fprintf('Open boundary ST forcing making time: %.2f minutes\n', toc / 60)

tic
% Write the temperature and salinity.
if strcmpi('HYCOM', {inputConf.obc_temp, inputConf.obc_salt})
    fprintf('Writing daliy open boundary S&T forcing file.\n')
    write_FVCOM_tsobc(fullfile(inputConf.outbase, inputConf.casename), ...
        Mobj.ts_times, ...
        size(Mobj.temperature, 2), ...
        Mobj.temperature, ...
        Mobj.salt, ...
        Mobj, ...
        'floattime', true,...
        'julian', true);
end

% Write the meanflow.
if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
    fprintf('Writing daliy open boundary meanflow file.\n')
    write_FVCOM_meanflow(Mobj, ...
        fullfile(inputConf.outbase,[inputConf.casename,'_mfobc.nc']), ...
        Mobj.velocity);
end

% plot temp and salt
% plot(datetime(Mobj.ts_times+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),squeeze(Mobj.temperature(1,1,:)));
% plot(datetime(Mobj.ts_times+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),squeeze(Mobj.salt(1,1,:)));

%%
%%%------------------------------------------------------------------------
%%%                     Meteorology and output
%%%------------------------------------------------------------------------
% Get the surface heating data.
if strcmpi(inputConf.doForcing, 'NCEP')
    
    % Use the OPeNDAP NCEP script to get the following parameters:
    %     - Potential evaporation rate (pevpr)[W/m^2]                   : for extracting land mask
    %     - u wind component (uwnd)[m/s]                                : for wind
    %     - v wind component (vwnd)[m/s]                                : for wind
    %     - Precipitation rate (prate)[Kg/m^2/s]                        : for precipitation
    %     - Downward Longwave Radiation Flux at Surface (dlwrf) [W/m^2] : for heat flux
    %     - Downward Solar Radiation Flux at Surface (dswrf) [W/m^2]    : for heat flux
    %     - Upward Longwave Radiation Flux at Surface (ulwrf) [W/m^2]   : for heat flux
    %     - Upward Solar Radiation Flux at Surface (uswrf) [W/m^2]      : for heat flux
    %     - Sea level pressure (pres) [Pa]                              : for air pressure
    %     - Air temperature at 2 m (air) [Kelvins]
    %     - Relative Humidity on Pressure Levels (rhum) [%]
    % The script also calculate the following parameters:
    %     - Momentum flux (tau)
    %     - Net solar radiation surface (nswrs = uswrf - dswrf)
    %     - Net longwave radiation surface (nlwrs = ulwrf - dlwrf)

    if develop_mode == 3
    	fprintf('Loading Model objet file...\n')
        load('varb/Mobj_03.mat');
        load('varb/forcing_interp.mat');
        fprintf('Done!\n');
    else
    	if develop_mode == 1
        	% The script converts the NCEP data from the OPeNDAP server from longitudes 
        	% in the 0 to 360 range to the latitudes in the -180 to 180 range. 
        	% It also subsets for the right region (defined by Mobj.lon and Mobj.lat).
        	% Uncomment the variables in get_NCEP_forcing as the varlist shows.
        	fprintf('Downloading NCEP forcing from OPeNDAP server database...\n')
            % modelTime=inputConf.forceMJD; varargin
        	forcing_ncep = get_NCEP_forcing(Mobj, inputConf.forceMJD, ...
        	    'varlist', {...
        	    'uwnd', 'vwnd',...
        	    'uswrf', 'ulwrf', 'dswrf', 'dlwrf',...
        	    'prate', 'pres','air', 'rhum','pevpr'},...
        	    'source', 'reanalysis2');
        	forcing_ncep.domain_cols = length(forcing_ncep.lon);
        	forcing_ncep.domain_rows = length(forcing_ncep.lat);
        	if isfield(forcing_ncep, 'rhum')||isfield(forcing_ncep, 'pres')
        	    forcing_ncep.domain_cols_alt = length(forcing_ncep.rhum.lon);
        	    forcing_ncep.domain_rows_alt = length(forcing_ncep.rhum.lat);
        	end
        	% Convert the small subdomain into cartesian coordinates. We need
        	% to do this twice because some of the NCEP data are on different
        	% grids (e.g. sea level pressure, relative humidity etc.).
        	tmpZone = regexpi(inputConf.utmZone,'\ ','split');
        	[tmpLon, tmpLat] = meshgrid(forcing_ncep.lon, forcing_ncep.lat);
        	[forcing_ncep.x, forcing_ncep.y] = wgs2utm(tmpLat(:), tmpLon(:), str2double(char(tmpZone{1}(1))), char(tmpZone{1}(2)));
        	if isfield(forcing_ncep, 'rhum')||isfield(forcing_ncep, 'pres')
        	    [tmpLon2, tmpLat2] = meshgrid(forcing_ncep.rhum.lon, forcing_ncep.rhum.lat);
        	    [forcing_ncep.xalt, forcing_ncep.yalt] = wgs2utm(tmpLat2(:), tmpLon2(:), str2double(char(tmpZone{1}(1))), char(tmpZone{1}(2)));
        	end
        	clear tmpLon tmpLat tmpLon2 tmpLat2 tmpZone
        	% Create arrays of the x and y positions.
        	forcing_ncep.x = reshape(forcing_ncep.x, forcing_ncep.domain_rows, forcing_ncep.domain_cols);
        	forcing_ncep.y = reshape(forcing_ncep.y, forcing_ncep.domain_rows, forcing_ncep.domain_cols);
        	if isfield(forcing_ncep, 'rhum')||isfield(forcing_ncep, 'pres')
        	    forcing_ncep.xalt = reshape(forcing_ncep.xalt, forcing_ncep.domain_rows_alt, forcing_ncep.domain_cols_alt);
        	    forcing_ncep.yalt = reshape(forcing_ncep.yalt, forcing_ncep.domain_rows_alt, forcing_ncep.domain_cols_alt);
        	end
        	[forcing_ncep.lon, forcing_ncep.lat] = meshgrid(forcing_ncep.lon, forcing_ncep.lat);
        	forcing_ncep = rmfield(forcing_ncep, {'domain_rows', 'domain_cols'});
        	if isfield(forcing_ncep, 'rhum')||isfield(forcing_ncep, 'pres')
        	    forcing_ncep = rmfield(forcing_ncep, {'domain_rows_alt', 'domain_cols_alt'});
        	end
        	fprintf('Saving NCEP forcing from OPeNDAP server database...')
        	save('forcing_ncep.mat','forcing_ncep','-v7.3','-nocompression');
        	fprintf('Done\n')
            % % Have a look at some data.
            %{
            for i=1:240%size(forcing_ncep.air.data, 3)
                figure1 = figure(1);
                ax = axes('Parent',figure1);
                clf
                %           1      2     3      4     ???5      6     7     8
                varb = {'uwnd','vwnd','air','rhum','prate','pres'};
                %value = sqrt(forcing_ncep.(varb{1}).data(:, :, i).^2 + forcing_ncep.(varb{2}).data(:, :, i).^2);
                value = forcing_ncep.(varb{3}).data(:,:,i);
                [X, Y] = meshgrid(forcing_ncep.lon, forcing_ncep.lat);
                s = pcolor(forcing_ncep.lon, forcing_ncep.lat, value);
                shading flat
                s.FaceColor = 'interp';
                axis('equal','tight')
                ylabel('Latitude (degree)','FontSize',12);
                xlabel('Longtitude (degree)','FontSize',12);
                pause(0.01);
            end
            clear ans ax figure1 i s value varb
            plot(squeeze(forcing_ncep.nswrs.data(2,4,1:36)));
            hold on;
            plot(squeeze(forcing_ncep.nlwrs.data(2,4,1:36)));
            hold on;
            plot(squeeze(forcing_ncep.shtfl.data(2,4,1:36)));
            hold on;
            plot(squeeze(forcing_ncep.lhtfl.data(2,4,1:36)));
            % The result of reanalysis2 shows the sensible and 
            % latent is totaly wrong.
            % reanalysis2 shtfl and lhtfl can not be used.
            %}
        elseif develop_mode == 2
    		fprintf('Loading NCEP forcing from the local database...');
        	load('forcing_ncep.mat');
        	fprintf('Done!\n');
    	end
        % Interpolate the data onto the FVCOM unstructured grid.
        interpfields = {'time', 'lon', 'lat', 'x', 'y',...
        	    'uwnd', 'vwnd',...
        	    'uswrf', 'ulwrf', 'dswrf', 'dlwrf',...
        	    'prate', 'pres', 'air', 'rhum','pevpr'};
        forcing_ncep_interp = grid2fvcom(Mobj, interpfields, forcing_ncep,...
            'add_elems', true);

        % save('varb/Mobj_03.mat','Mobj','-v7.3','-nocompression');
        % save('forcing_interp.mat','forcing_interp','-v7.3','-nocompression');
    end
    % elseif strcmpi(inputConf.doForcing, 'NCEP-CALCULATED')
    % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.uwnd.node(1,:));
    % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.vwnd.node(1,:));
    % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.prate.node(1,:));
    % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.evap.node(1,:));
    % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.nswrs.node(1,:));
    % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.dlwrf.node(1,:));
    % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.air.node(1,:));
    % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.pres.node(1,:));
    % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.rhum.node(1,:));
    % write_FVCOM_forcing_calculated(Mobj, ...
    %     fullfile(inputConf.outbase,inputConf.casename),...
    %     forcing_interp_calculated,...
    %     [inputConf.doForcing, 'atmospheric forcing data'],...
    %     inputConf.FVCOM_version, ...
    %     'floattime', true,...
    %     'julian', true);
    % fileprefix=fullfile(inputConf.outbase,inputConf.casename);
    % data=forcing_interp_calculated;
    % infos=[inputConf.doForcing, 'atmospheric forcing data'];
    % fver=inputConf.FVCOM_version;

elseif strcmpi(inputConf.doForcing, 'GWO')
    
    % Get the surface heating data.
    % Use the GWO data to get the following parameters:
    % - u wind component (uwnd)[m/s]                       : for wind
    % - v wind component (vwnd)[m/s]                       : for wind
    % - Precipitation rate (prate)[m/s]                    : for precipitation
    % - Solar Radiation Flux (dsw) [W/m^2]                 : for heat flux
    % - Sea level pressure (pres) [Pa]                     : for air pressure
    % - Air temperature at ? m (air) [deg C]               : for temperature
    % - Relative Humidity on Pressure Levels (rhum) [%]    : for relative humidity
    % - Cloud coverage (cloud)

    % Generate the following variables for FVCOM net sea surface heat flux and sea
    % surface shortwave radiation. For this purpose, heat_calculated and custom
    % created dlw_calculated options are actived.
    %  &NML_SURFACE_FORCING
    %  WIND_ON                         = T,
    %  WIND_TYPE                       = 'speed',
    %  WIND_FILE                       = 'tokyobay_hfx.nc', 
    %  &NML_HEATING_CALCULATED
    %  HEATING_CALCULATE_ON            = T,
    %  HEATING_CALCULATE_TYPE          = 'flux',
    %  HEATING_CALCULATE_FILE          = 'tokyobay_hfx.nc',
    %  COARE_VERSION                   = 'BULKALGORITHM',
    %  The air pressure is not turnning on.
    %  This is used for non calculated heat flux and hurrican model.
    
    if develop_mode == 3
    	fprintf('Loading Model objet file...\n')
        load('varb/Mobj_03.mat');
        % load('forcing_gwo.mat');
        load('forcing_gwo_interp.mat');
        fprintf('Done!\n');
    else
        Mobj.gwo.time = inputConf.forceMJD';
        if develop_mode == 1
            forcing_gwo = get_GWO_forcing(Mobj.gwo.time);
            UTMzone = regexpi(inputConf.utmZone,'\ ','split');
            for s = 1:length(forcing_gwo.lon)
                [forcing_gwo.x(s),forcing_gwo.y(s),~,~] = wgs2utm(...
                    forcing_gwo.lat(s),forcing_gwo.lon(s),...
                    str2double(char(UTMzone{1}(1))),char(UTMzone{1}(2)));
            end
            [forcing_gwo.x,forcing_gwo.y] = meshgrid(forcing_gwo.x,forcing_gwo.y);
            clear s UTMzone
            save('forcing_gwo.mat','forcing_gwo','-v7.3','-nocompression');
            save('varb/Mobj_03.mat','Mobj','-v7.3','-nocompression');
        elseif develop_mode == 2
    		fprintf('Loading GWO forcing from the local database...');
        	load('forcing_gwo.mat');
        	fprintf('Done!\n');
        end
        % air
        interpfields = {'time', 'lon', 'lat', 'x', 'y',...
            'air'};
        % vars = interpfields; data = forcing;
        forcing_gwo_inter_air = grid2fvcom(Mobj, interpfields, forcing_gwo,...
            'add_elems', false);
        save('forcing_gwo_inter_air.mat','forcing_gwo_inter_air','-v7.3','-nocompression');
        clear forcing_gwo_inter_air;
        
        % cld
        interpfields = {'time', 'lon', 'lat', 'x', 'y',...
            'cld'};
        % vars = interpfields; data = forcing;
        forcing_gwo_inter_cld = grid2fvcom(Mobj, interpfields, forcing_gwo,...
            'add_elems', false);
        save('forcing_gwo_inter_cld.mat','forcing_gwo_inter_cld','-v7.3','-nocompression');
        clear forcing_gwo_inter_cld;
        
        % dsw
        interpfields = {'time', 'lon', 'lat', 'x', 'y',...
            'dsw'};
        % vars = interpfields; data = forcing;
        forcing_gwo_inter_dsw = grid2fvcom(Mobj, interpfields, forcing_gwo,...
            'add_elems', false);
        save('forcing_gwo_inter_dsw.mat','forcing_gwo_inter_dsw','-v7.3','-nocompression');
        clear forcing_gwo_inter_dsw;
        
        % hum
        interpfields = {'time', 'lon', 'lat', 'x', 'y',...
            'hum'};
        % vars = interpfields; data = forcing;
        forcing_gwo_inter_hum = grid2fvcom(Mobj, interpfields, forcing_gwo,...
            'add_elems', false);
        save('forcing_gwo_inter_hum.mat','forcing_gwo_inter_hum','-v7.3','-nocompression');
        clear forcing_gwo_inter_hum;
        
        % prs
        interpfields = {'time', 'lon', 'lat', 'x', 'y',...
            'prs'};
        % vars = interpfields; data = forcing;
        forcing_gwo_inter_prs = grid2fvcom(Mobj, interpfields, forcing_gwo,...
            'add_elems', false);
        save('forcing_gwo_inter_prs.mat','forcing_gwo_inter_prs','-v7.3','-nocompression');
        clear forcing_gwo_inter_prs;
        
        % wnd
        interpfields = {'time', 'lon', 'lat', 'x', 'y',...
            'uwnd','vwnd'};
        % vars = interpfields; data = forcing_gwo;
        forcing_gwo_inter_wnd = grid2fvcom(Mobj, interpfields, forcing_gwo,...
            'add_elems', true);
        save('forcing_gwo_inter_wnd.mat','forcing_gwo_inter_wnd','-v7.3','-nocompression');
        clear forcing_gwo_inter_wnd;
        
        % rin
        interpfields = {'time', 'lon', 'lat', 'x', 'y',...
            'rin'};
        % vars = interpfields; data = forcing;
        forcing_gwo_inter_rin = grid2fvcom(Mobj, interpfields, forcing_gwo,...
            'add_elems', false);
        save('forcing_gwo_inter_rin.mat','forcing_gwo_inter_rin','-v7.3','-nocompression');
        clear forcing_gwo_inter_rin;
        clear forcing_gwo;
        clear interpfields;
        
        load('forcing_gwo_inter_air.mat'); forcing_gwo_interp.air = forcing_gwo_inter_air.air; clear forcing_gwo_inter_air;
        load('forcing_gwo_inter_cld.mat'); forcing_gwo_interp.cld = forcing_gwo_inter_cld.cld; clear forcing_gwo_inter_cld;
        load('forcing_gwo_inter_dsw.mat'); forcing_gwo_interp.dsw = forcing_gwo_inter_dsw.dsw; clear forcing_gwo_inter_dsw;
        load('forcing_gwo_inter_hum.mat'); forcing_gwo_interp.hum = forcing_gwo_inter_hum.hum; clear forcing_gwo_inter_hum;
        load('forcing_gwo_inter_prs.mat'); forcing_gwo_interp.prs = forcing_gwo_inter_prs.prs; clear forcing_gwo_inter_prs;
        load('forcing_gwo_inter_rin.mat'); forcing_gwo_interp.rin = forcing_gwo_inter_rin.rin; clear forcing_gwo_inter_rin;
        load('forcing_gwo_inter_wnd.mat'); 
        forcing_gwo_interp.uwnd = forcing_gwo_inter_wnd.uwnd; 
        forcing_gwo_interp.vwnd = forcing_gwo_inter_wnd.vwnd; 
        forcing_gwo_interp.time = forcing_gwo_inter_wnd.time;
        forcing_gwo_interp.lon  = forcing_gwo_inter_wnd.lon; 
        forcing_gwo_interp.lat  = forcing_gwo_inter_wnd.lat; 
        forcing_gwo_interp.x    = forcing_gwo_inter_wnd.x; 
        forcing_gwo_interp.y    = forcing_gwo_inter_wnd.y; 
        clear forcing_gwo_inter_wnd;
        
        save('forcing_gwo_interp.mat','forcing_gwo_interp','-v7.3','-nocompression');
    end

    write_FVCOM_gwo_forcing(Mobj, ...
        forcing_gwo_interp, ...
        fullfile(inputConf.outbase,[inputConf.casename]),...
        [inputConf.doForcing, ' atmospheric forcing data'],...
        inputConf.FVCOM_version, ...
        'floattime', true,...
        'julian', true);
    fprintf('Writing meteorological forcing file...done!\n');
    clear forcing_gwo forcing_gwo_interp;
    
% Have a look at some data.

for i=1:240%size(forcing.air.data, 3)
    figure1 = figure(1);
    ax = axes('Parent',figure1);
    clf
    %           1      2     3     4     5     6     7     8
    varb = {'uwnd','vwnd','air','hum','rin','prs','cld','dsw'};
    %value = sqrt(forcing_gwo.(varb{1}).data(:, :, i).^2 + forcing_gwo.(varb{2}).data(:, :, i).^2);
    value = forcing_gwo.(varb{8}).data(:,:,i);
    s = pcolor(forcing_gwo.lon, forcing_gwo.lat, value');
    shading flat
    s.FaceColor = 'interp';
    axis('equal','tight')
    ylabel('Latitude (degree)','FontSize',12);
    xlabel('Longtitude (degree)','FontSize',12);
    pause(0.01);
end
clear ans ax figure1 i s value varb

end
end

%%
%%%------------------------------------------------------------------------
%%%                     River discharge and output
%%%------------------------------------------------------------------------
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
        save('varb/Mobj_04.mat','Mobj','-v7.3','-nocompression');
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
diary off;


%load('/Users/yulong/GitHub/water/data/monitoring_post/case_2014_2017/layers/water.mat');

%plot(datetime(water.sf.DATETIME(1:7932),'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),water.sf.TEMPdeg(1:7932));
%plot(datetime(water.md.DATETIME(1:7932),'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),water.md.TEMPdeg(1:7932));
%plot(datetime(water.bt.DATETIME(1:7747),'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),water.bt.TEMPdeg(1:7747));

%plot(datetime(water.sf.DATETIME(1:7932),'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),water.sf.SALpsu(1:7932));
%plot(datetime(water.md.DATETIME(1:7932),'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),water.md.SALpsu(1:7932));
%plot(datetime(water.bt.DATETIME(1:7747),'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),water.bt.SALpsu(1:7747));
