function Mobj= get_COSTUM_river_variable(Mobj, fileprefix, name, infos, varargin)

%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\Yulong WANG\Documents\GitHub\river\data_q\newdata_2014_2017.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

if (nargin < 2 || nargin > 3) && isempty(varargin)
    error('Incorrect number of arguments')
end

% By Yulong Wang on 2019/09/09

%%
for i =1:length(fileprefix)
    %% Initialize variables.
    filename = fileprefix{i};
    delimiter = ',';
    startRow = 2;
    riverNum = length(infos);
    % Default to string times as FVCOM looks for these first.
    time = false;
    for vv = 1:2:length(varargin)
        switch varargin{vv}
            case 'time'
                time = true;
        end
    end

    %% Format for each line of text:
    %   column1: datetimes (%{MM/dd/yyyy HH:mm}D)
    %	column2: double (%f)
    %   column3: double (%f)
    %	column4: double (%f)
    %   column5: double (%f)
    %	column6: double (%f)
    %   column7: double (%f)
    %	column8: double (%f)
    % For more information, see the TEXTSCAN documentation.
    % formatSpec = '%{MM/dd/yyyy HH:mm}D%f%f%f%f%f%f%f%[^\n\r]';
    formatSpec = ['%{MM/dd/yyyy HH:mm}D',repmat('%f',1,riverNum),'%[^\n\r]'];
    % formatSpec = ['%{yyyy/mm/dd HH:mm}D',repmat('%f',1,riverNum),'%[^\n\r]'];

    %% Open the text file.
    fileID = fopen(filename,'r');

    %% Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

    %% Close the text file.
    fclose(fileID);

    %% Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post
    % processing code is included. To generate code which works for
    % unimportable data, select unimportable cells in a file and regenerate the
    % script.

    %% Create output variable
    data = table(dataArray{1:end-1}, 'VariableNames', ['Time',infos]);

    % For code requiring serial dates (datenum) instead of datetime, uncomment
    % the following line(s) below to return the imported dates as datenum(s).
    % for example 01/01/2014 00:00 > 735600

    data.Time=datenum(data.Time);

    % Transformation from table to array.
    array = table2array(data);

    Mobj.river.(name{i}) = array(:,2:end);

end

if time
    Mobj.river.time = array(:,1);
    Mobj.river.Time = datestr(Mobj.river.time,31);
    Mobj.river.timeMJD = array(:,1) - 678942;
    fprintf('\nWARNING: River time begin from %s\n', Mobj.river.Time(1,1:end));
end

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%%
return