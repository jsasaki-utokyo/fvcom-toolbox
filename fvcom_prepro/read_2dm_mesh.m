function [Mobj] = read_2dm_mesh(varargin) 

% Read .grid mesh files into Matlab mesh object  
%
% [Mobj] = function read_2dm_mesh(varargin)
%
% DESCRIPTION:
%    Read NOCL .grid file 
%    Store in a matlab mesh object 
%
% INPUT [keyword pairs]:  
%   'grid'                  = NOCL .grid file (e.g. UK_grid.2dm)
%   [optional] 'coordinate' = coordinate system for output data [spherical; cartesian (default)]
%   [optional] 'input_coord' = coordinate system for input data [spherical; cartesian (default)]
%   [optional] 'project'    = generate (x,y) coordinates if input is (lon,lat) 
%                             generate (lon,lat) coordinates if input is (x,y)
%                            ['true' ; 'false' (default)], see myproject.m
%   [optional] 'zone'    = specify UTM zone for projection
%   [optional] 'addCoriolis' = calculate Coriolis param (f), requires [lon,lat]
%
% OUTPUT:
%    Mobj = matlab structure containing mesh data
%
% EXAMPLE USAGE
%    Mobj = read_grid_mesh('2dm','UK_grid.2dm','coordinate','spherical','addCoriolis','true')
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Thurston (National Oceanography Centre Liverpool)
%
% Revision history (KJT)
%   2012-11-13 Adapted 'read_sms_mesh.m' to read NOCL .grid files
%   2012-11-19 Added input/output coordinate functionality
%   
%==============================================================================

subname = 'read_grid_mesh';
global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

have_2dm = false;
have_bath = false;
have_lonlat = false;
have_xy = false;
userproject = false;
haveUTM = false;
addCoriolis = false;

%------------------------------------------------------------------------------
% Create a blank mesh object
%------------------------------------------------------------------------------
Mobj = make_blank_mesh();

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------

if(mod(length(varargin),2) ~= 0)
	error(['incorrect usage of ',subname,', use keyword pairs'])
end

for i=1:2:length(varargin)-1
	keyword  = lower(varargin{i});
	if( ~ischar(keyword) )
		error(['incorrect usage of ',subname,', check keywords'])
    end
	
	switch(keyword(1:3))
	
    case '2dm'
		sms2dm_fname = varargin{i+1};
		have_2dm = true;
	case 'coo'
		coord = varargin{i+1};
		if(coord(1:3)=='sph')
			coordinate = 'spherical';
		else
			coordinate = 'cartesian';
        end
    case 'in_'
		in_coord = varargin{i+1};
		if(in_coord(1:3)=='sph')
			in_coordinate = 'spherical';
            have_lonlat = true;
		else
			in_coordinate = 'cartesian';
            have_xy = true;
        end
	case 'pro'
		val = varargin{i+1};
		if(val)
			userproject = true;
		else
			userproject = false;
        end
    case 'zon'
        fullzone = varargin{i+1};
        UTMzone = regexpi(fullzone,'\ ','split');
        UTMzone=str2double(char(UTMzone{1}(1)));
        haveUTM = true;
    case 'add'
		val = varargin{i+1};
		if( val )
			addCoriolis = true;
        else
			addCoriolis = false;
        end
	otherwise
		error(['Can''t understand property:' char(varargin{i+1})]);
    end %switch(keyword)
end
		
%------------------------------------------------------------------------------
% Read the mesh from the .2dm file
%------------------------------------------------------------------------------

fid = fopen(sms2dm_fname,'r');
if(fid  < 0)
	error(['file: ' sms2dm_fname ' does not exist']);
end

% Count number of elements and vertices
if(ftbverbose)
  fprintf(['reading from: ' sms2dm_fname '\n'])
end

nElems = 0;
nVerts = 0;
nStrings = 0;
nHeader = 0;

StillReading = true;
while StillReading
    lin = fgetl(fid);
    if lin ~= -1 % EOF is -1
        switch lin(1:2)
            case 'E3'
                nElems = nElems + 1;
            case 'ND'
                nVerts = nVerts + 1;
            case 'NS'
                nStrings = nStrings + 1;
            case {'ME', 'NU'}
                nHeader = nHeader + 1;
            case 'E4'
                error('Quadrilateral elements are unsupported in FVCOM')
            otherwise
                StillReading = false;
        end
    else
        % Got to EOF
        StillReading = false;
    end
end
fclose(fid);

fid = fopen(sms2dm_fname, 'rt');

if(ftbverbose)
  fprintf('nVerts: %d\n',nVerts); 
  fprintf('nElems: %d\n',nElems); 
  fprintf('reading in connectivity and grid points\n')
end

% allocate memory to hold mesh and connectivity
tri = zeros(nElems,3);
lon = zeros(nVerts,1);
lat = zeros(nVerts,1);
ts  = zeros(nVerts,1);

% Skip the header
% C = textscan(fid, '%s', nHeader + 1,'Headerlines',nHeader); 
% why does it skip 3 lines when only two are in the header? Didn't use to throw errors though...
% if Meshname is more than one word the above would fail... Using
% headerlines as below is more safe proof.
% Read the triangulation table
C = textscan(fid, '%s %d %d %d %d %d', nElems,'Headerlines',nHeader);
tri(:, 1) = C{3};
tri(:, 2) = C{4};
tri(:, 3) = C{5};

% Read in the nodes and interpolated depths
C = textscan(fid, '%s %d %f %f %f ', nVerts);
x = C{3};
y = C{4};
h = C{5};


% Check we don't have any NaNs anywhere
if max(isnan(x)) == 1
    error('%d NaNs in the x data', sum(isnan(x)))
end
if max(isnan(y)) == 1
    error('%d NaNs in the y data', sum(isnan(x)))
end
if max(isnan(h)) == 1
    error('%d NaNs in the h data', sum(isnan(x)))
end
if max(isnan(tri(:))) == 1
    error('%d NaNs in the h data', sum(isnan(tri(:))))
end

% Build array of all the nodes in the nodestrings.
C = textscan(fid, '%s %d %d %d %d %d %d %d %d %d %d', nStrings);
allNodes = cell2mat(C(2:end))';
allNodes(allNodes == 0) = [];

% % Add a new field to Mobj with all the nodestrings as a cell array.
nodeStrings = find(allNodes < 0);
read_obc_nodes = cell(1, length(nodeStrings));
for nString = 1:sum(allNodes(:) < 0)
    if nString == 1
        read_obc_nodes{nString} = abs(allNodes(1:nodeStrings(nString)));
    else
        read_obc_nodes{nString} = ...
            abs(allNodes(nodeStrings(nString - 1) + 1: ...
            nodeStrings(nString)));
    end
    % Check for closed nodestrings (which we don't really want).
    if read_obc_nodes{nString}(1) == read_obc_nodes{nString}(end)
        % Drop the end one.
        warning('Closed node string found. Opening it.')
        read_obc_nodes{nString} = read_obc_nodes{nString}(1:end - 1);
    end
end

if nStrings > 0
    have_strings = true;
end

if have_lonlat
    lon = x;
    lat = y;
    % Just do a double check on the coordinates to make sure we don't
    % actually have cartesian
    if max(lon) > 360
        warning('You''ve specified spherical coordinates, but your upper longitude value exceeds 360 degrees. Are you sure you have spherical data?')
    end
elseif have_xy
    have_xy = true;
else
    warning('Unrecognised coordinate system (%s). Valid values are ''spherical'' and ''cartesian''.', coordinate)
end

fclose(fid);

% Make sure we have bathymetry
if sum(h)==0
    have_bath=false;
else
    have_bath=true;
end

% Make sure we have positive depths
if sum(h>0) < sum(h<0)
    h = -h;
end

%------------------------------------------------------------------------------
% Project if desired by user
%------------------------------------------------------------------------------

if(userproject)
    if (in_coordinate(1:3)=='car')
        fprintf('inverse projecting to get (lon,lat)\n')
        utmZones=cellfun(@(x) repmat(x,length(x),1),fullzone,'uni',false);
        [lon,lat] = utm2deg(x,y,utmZones{1});
        Mobj.have_lonlat = true;
    elseif (in_coordinate(1:3)=='sph')
        fprintf('forward projecting to get (x,y)\n')
        [x,y] = wgs2utm(lat,lon,UTMzone,'N');
        have_xy = true;
    end
end

%------------------------------------------------------------------------------
% Transfer to Mesh structure
%------------------------------------------------------------------------------

Mobj.nVerts  = nVerts;
Mobj.nElems  = nElems;
Mobj.nativeCoords = coordinate;

if have_lonlat==true
	Mobj.have_lonlat    = have_lonlat;
end
if have_xy==true
	Mobj.have_xy        = have_xy;
end

Mobj.have_bath      = have_bath;

Mobj.read_obc_nodes = read_obc_nodes;

Mobj.x            = x;
Mobj.y            = y;
Mobj.ts           = ts;
Mobj.lon          = lon;
Mobj.lat          = lat;
Mobj.h            = h;
Mobj.tri          = tri;

if addCoriolis==true
    Mobj.have_cor = true;
    Mobj = add_coriolis(Mobj,'uselatitude');
end

if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end


