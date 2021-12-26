function header=read_header(entry,words)
% sets the header values according to what was read in the .ctl-file.
% cf. GrADS documentation : http://grads.iges.org/grads/gadoc/descriptorfile.html

global l_quiet

header=evalin('caller','header');
fid_ctl=header.FID{1};

grid_fields={'num','type','vec','rev'};

switch upper(entry)
 case 'DSET'                    % Name of the bin. data file.
  param=deal(words{1});
  if param(1)=='^' & isfield(header,'DIR')
    param(1)=[];
    param=[header.DIR,param];
  elseif param(1)=='^'
    param(1)=[];
  end
  if exist(param,'file')
    header.DATANAME=param; 
  else
    error(['Non-existent file ',param])
  end
  %
 case 'DTYPE'                  % Gridded or station data-file (not used later)
  param=lower(words{1});
  switch param
   case 'grib'
    header.DTYPE=1;
   case 'station'
    header.DTYPE=0;
  end
  %
 case 'STNMAP'                  % Name of the map file for station data (set in INDEX)
  if header.DTYPE
    warning([header.FILENAME,' cannot specify STNMAP=',param,' for a gridded data set.'])
  end
  param=deal(words{1});
  if param(1)=='^'
    param(1)=[];
  end
  if exist(param,'file')
    header.INDEX=param; 
  else
    error(['Non-existent file ',param])
  end
  %
 case 'INDEX'                  % Name of the GRIB map file (not used later)
  [header.INDEX]=strcat(words{:});
  %
 case 'TITLE'                  % Data title (not used later)
  line={words{:}};
  space=cell(size(line));
  [space{:}]=deal(char(32));
  line=cat(1,line,space);
  line=reshape(line,1,prod(size(line)));   % Rearranging the title to keep spaces between words.
  header.TITLE=deblank(cat(2,line{:}));
  %
 case 'UNDEF'                  % Missing or undefined data value
  param=words{1};
  if isreal(str2num(param))
    header(1).UNDEF=str2num(param);
  else
    error(['Wrong missing value format : ',param])
  end
  %
 case 'OPTIONS'                % Possible options (not all implemented...)
  for param=words
  param=param{1};
  if ~ischar(param), error(['Wrong option format : ',param]), end
  switch lower(param)
   case 'big_endian'
    header.BINTYPE='ieee-be';
   case 'cray_32bit_ieee'
    header.BINTYPE='cray';
   case 'little_endian'
    header.BINTYPE='ieee-le';
   case 'binprecision'
    % Strange syntax to eliminate digits in BINPRECISION
    param=words{2}(find(isnan(str2double(cellstr(upper(words{2})')))));
    if any(strmatch(param,...
		    {'uchar','schar','int','uint','float','single','double',...
		     'real*','integer*','bit','ubit','char','short','ushort', ...
		     'long','ulong'}))
      header.BINPRECISION=words{2};
    end
   case 'xrev'
    header.XDEF.rev=1;
   case 'yrev'
    header.YDEF.rev=1;
   case 'zrev'
    header.ZDEF.rev=1;
   case 'sequential' 
    % Structure of the binary file: ([4 bit][XYHEADER bit][XDEF.num
    % x YDEF.num x BINPRECISION bit][4 bit]) x NVAR
    header.SEQ=1;
   otherwise
    fprintf('Option %s not yet supported.\n',param);
  end
  end % for param=words
  %
% $$$   Matlab     C/Fortran       Description  
% $$$   'uchar'   'unsigned char'  unsigned character,  8 bits.
% $$$   'schar'   'signed char'    signed character,  8 bits.
% $$$   'int8'    'integer*1'      integer, 8 bits.
% $$$   'int16'   'integer*2'      integer, 16 bits.
% $$$   'int32'   'integer*4'      integer, 32 bits.
% $$$   'int64'   'integer*8'      integer, 64 bits.
% $$$   'uint8'   'integer*1'      unsigned integer, 8 bits.
% $$$   'uint16'  'integer*2'      unsigned integer, 16 bits.
% $$$   'uint32'  'integer*4'      unsigned integer, 32 bits.
% $$$   'uint64'  'integer*8'      unsigned integer, 64 bits.
% $$$   'single'  'real*4'         floating point, 32 bits.
% $$$   'float32' 'real*4'         floating point, 32 bits.
% $$$   'double'  'real*8'         floating point, 64 bits.
% $$$   'float64' 'real*8'         floating point, 64 bits.

 case 'XDEF'                   % Description of the X-axis (longitude usually)
  % Checking the number of x-grid coordinates:
  if isreal(str2num(words{1}))
    words{1}=str2num(words{1});
  elseif str2num(words{1})>1
  else
    error(['Wrong number of x-grid arguments : ', ...
					num2str(words{1})])
  end
  
  switch upper(words{2})
   case 'LINEAR'
    vec(1)=str2num(words{3});
    vec(2)=str2num(words{4});
    vec2=vec(1) + vec(2) * [0:(words{1}-1)];
    words{3}=vec2;				% KS, 24.11.2003: coord vector in all cases. 
   case 'LEVELS'
% $$$                 try % all on the same line
% $$$                     for i=1:words{1}
% $$$                         vec(i)=str2num(words{i+2});
% $$$                     end
% $$$ 		    words{3}=vec;
% $$$                 catch % all levels written column-wise
% $$$                     for i=1:words{1}
% $$$                         % vec(i)=fscanf(fid_ctl,'%i',1);
% $$$                         vec(i)=str2num(fgetl(fid_ctl));
% $$$                     end
% $$$ 		    words{3}=vec;
% $$$                 end
%   vec=str2num([words{3:end}]');
%   vec=str2num(cell2mat({words{3:end}}'))';
    vec=str2num(strvcat(words{3:end}));
    
    while length(vec) < words{1}
      point=ftell(fid_ctl);
      read_line=str2num(fgetl(fid_ctl));
      if isempty(read_line)
	fseek(fid_ctl,point,'bof');
	break
      else
	vec=cat(2,vec,read_line);
      end
    end
    if length(vec)==words{1}
      words{3}=vec;
    else
      error(['Unsuitable number of X-levels : ',num2str(length(vec))])
    end
   otherwise
    error(['Unsuitable x-mapping method : ',words{2}])
  end
  
  %	if header(1).XDEF.rev, words{3}=fliplr(words{3}); end
  [header(1).XDEF]=cell2struct({words{1},words{2},words{3},header(1).XDEF.rev},grid_fields,2);
  
  %
 case 'YDEF'                   % Description of the Y-axis (latitude usually)
   % Checking the number of y-grid coordinates:
   if isreal(str2num(words{1}))
     words{1}=str2num(words{1});
   elseif str2num(words{1})>1
   else
     error(['Wrong number of y-grid arguments : ', ...
	    num2str(words{1})])
   end
   
   switch upper(words{2})
    case 'LINEAR'
     vec(1)=str2num(words{3});
     vec(2)=str2num(words{4});
     vec2=vec(1) + vec(2) * [0:(words{1}-1)];
     words{3}=vec2;
    case 'LEVELS'
% $$$                 try % all on the same line
% $$$                     for i=1:words{1}
% $$$                         vec(i)=str2num(words{i+2});
% $$$                     end
% $$$ 		    words{3}=vec;
% $$$                 catch % all levels written column-wise
% $$$                     for i=1:words{1}
% $$$                         % vec(i)=fscanf(fid_ctl,'%i',1);
% $$$                         vec(i)=str2num(fgetl(fid_ctl));
% $$$                     end
% $$$ 		    words{3}=vec;
% $$$                 end
%    vec=str2num([words{3:end}]');
%    vec=str2num(cell2mat({words{3:end}}'))';
     vec=str2num(strvcat(words{3:end}));
     while length(vec) < words{1}
       point=ftell(fid_ctl);
       read_line=str2num(fgetl(fid_ctl));
       if isempty(read_line)
	 fseek(fid_ctl,point,'bof');
	 break
       else
	 vec=cat(2,vec,read_line);
       end
     end
     if length(vec)==words{1}
       words{3}=vec;
     else
       error(['Unsuitable number of Y-levels : ',num2str(length(vec))])
     end
    otherwise
     error(['Unsuitable y-mapping method : ',words{2}])
   end
   
%	if header.YDEF.rev, words{3}=fliplr(words{3}); end
   [header(1).YDEF]=cell2struct({words{1},words{2},words{3},header(1).YDEF.rev},grid_fields,2);
   %
 case 'ZDEF'                   % Description of the Z-axis
    % Checking the number of z-grid coordinates:
    if isreal(str2num(words{1}))
      words{1}=str2num(words{1});
    elseif str2num(words{1})>1
    else
      error(['Wrong number of z-grid arguments : ', ...
	     num2str(words{1})])
    end
    
    switch upper(words{2})
     case 'LINEAR'
      vec(1)=str2num(words{3});
      vec(2)=str2num(words{4});
      vec2=vec(1) + vec(2) * [0:(words{1}-1)];
      words{3}=vec2;
     case 'LEVELS'
% $$$                 try % all on the same line
% $$$                     for i=1:str2num(words{1})
% $$$                         vec(i)=str2num(words{i+2});
% $$$                     end
% $$$ 		    words{3}=vec;
% $$$                 catch % all levels written column-wise
% $$$                     for i=1:words{1}
% $$$                         % vec(i)=fscanf(fid_ctl,'%i',1);
% $$$                         vec(i)=str2num(fgetl(fid_ctl));
% $$$                     end
% $$$ 		    words{3}=vec;
% $$$                 end
%     vec=str2num([words{3:end}]');
%     vec=str2num(cell2mat({words{3:end}}'))';
      vec=str2num(strvcat(words{3:end}));
      while length(vec) < words{1}
	point=ftell(fid_ctl);
	read_line=str2num(fgetl(fid_ctl));
	if isempty(read_line)
	  fseek(fid_ctl,point,'bof');
	  break
	else
	  vec=cat(2,vec,read_line);
	end
      end
      if length(vec)==words{1}
	words{3}=vec;
      else
	error(['Unsuitable number of Z-levels : ',num2str(length(vec))])
      end
     otherwise
      error(['Unsuitable z-mapping method : ',words{2}])
    end
    
    %	if header.ZDEF.rev, words{3}=fliplr(words{3}); end
    [header(1).ZDEF]=cell2struct({words{1},words{2},words{3},header(1).ZDEF.rev},grid_fields,2);
    %
 case 'TDEF'                   % Description of the T-axis
  % Checking the number of t-grid coordinates:
  if isreal(str2num(words{1}))
    words{1}=str2num(words{1});
  elseif str2num(words{1})>1
  else
    error(['Wrong number of t-scale argument : ', ...
	   num2str(words{1})])
  end
  
  if any(upper(words{2})~='LINEAR')
    error(['Unsuitable t-mapping method : ',words{2}])
  end
  
  param=words{3};
  if (param(3)~=':' | param(6)~='Z' | ...
      ~isreal(str2num(param([1:2 4:5 7:8 12:length(param)]))))
    error(['Unsuitable date format : ',param])
  else
    vec(1)=str2num(param(1:2));          % hour
    vec(2)=str2num(param(4:5));          % min
    vec(3)=str2num(param(7:8));          % day
    switch lower(param(9:11))            % month
     case 'jan'
      vec(4)=1;
     case 'feb'
      vec(4)=2;
     case 'mar'
      vec(4)=3;
     case 'apr'
      vec(4)=4;
     case 'may'
      vec(4)=5;
     case 'jun'
      vec(4)=6;
     case 'jul'
      vec(4)=7;
     case 'aug'
      vec(4)=8;
     case 'sep'
      vec(4)=9;
     case 'oct'
      vec(4)=10;
     case 'nov'
      vec(4)=11;
     case 'dec'
      vec(4)=12;
     otherwise
      error(['Wrong month format : ',param(9:11)])
    end
    
    vec(5)=str2num(param(12:length(param)));       % year
    if vec(5)<10
      vec(5)=vec(5)+2000;
    elseif vec(5)<100
      vec(5)=vec(5)+1900;
    end
    
    % The time-axis is written in the matlab-time format
    % (cf. datestr, datevec...):
    init_time=datenum(vec(5),vec(4),vec(3),vec(1),vec(2),0);
    
    param=words{4};
    if isreal(str2num(param(1:end-2)))
      time_increment=str2num(param(1:end-2));
      param(1:end-2)=[];
    else
      error(['Wrong time increment argument : ',param])
    end
    
    switch lower(param)
     case 'mn'
      time_increment=time_increment*datenum(0,0,0,0,1,0);
     case 'hr'
      time_increment=time_increment*datenum(0,0,0,1,0,0);
     case 'dy'       % NB: datenum(0,0,1,0,0,0) = 1.
      time_increment=time_increment;
     case 'mo'       % standardised month
      time_increment=time_increment*365.25/12;
     case 'yr'       % standardised year
      time_increment=time_increment*365.25;
     otherwise
      error(['Wrong time increment argument : ',param])
    end
    
    % Writing the final time axis:
    clear vec
    vec(1)=init_time;
    vec(2)=time_increment;
    
    [header(1).TDEF]=cell2struct({words{1},words{2},vec,0},grid_fields,2);
  end
  %
 case 'VARS'                   % Number and size of variables
  if isreal(str2num(words{1}))
    n_var = str2num(words{1});
  else
    error(['Wrong variable number argument : ',words{1}])
  end
  
  % Reading the variables : 
  
  var_format={'id','name','levs','units','descr'};
  varl={};
  
  for i=1:n_var
    line_var=fgetl(fid_ctl);
    n_blank=length(line_var);
    line_var=strrep(line_var,'  ',' ');
    while length(line_var) < n_blank
      n_blank=length(line_var);
      line_var=strrep(line_var,'  ',' ');
    end
    while isspace(line_var(1)), line_var(1)=[]; end
    var_loc={i,'',NaN,'',''};
    var_hole=min(find(isspace(line_var)));
    var_loc{2}=line_var(1:(var_hole-1));
    line_var(1:var_hole)=[];
    
    var_hole=min(find(isspace(line_var)));
    var_loc{3}=str2num(line_var(1:(var_hole-1)));
  % cf GrADS doc:  If levs is 0, the variable does not correspond
  % to any vertical level
    var_loc{3}=max([var_loc{3},1]);
    line_var(1:var_hole)=[];
    
    var_hole=min(find(isspace(line_var)));
    var_loc{4}=line_var(1:(var_hole-1));
    line_var(1:var_hole)=[];
    
    var_loc{5}=line_var;
    
    varl=cat(1,varl,var_loc);
  end
  
  vars=cell2struct(varl,var_format,2);
  
  if ~l_quiet
    disp(['Is it the end ? ',fgets(fid_ctl)]);
  else
    endvars=deblank(fgets(fid_ctl));
    if ~strcmpi((endvars),'ENDVARS')
      error(['Wrong end of file : ',endvars])
    end
  end
  
  header(1).VARS=vars;
  assignin('caller','vars',vars);
  %
 case 'FILEHEADER'              %  n-byte header of binary data, not to be read
  [header.FILEHEADER]=str2num(words{1});
  %
 case 'THEADER'                 % n-byte header for each T-block, not to be read
  [header.THEADER]=str2num(words{1});
  %
 case 'XYHEADER'                % n-byte header for each XY-block, not to be read
  [header.XYHEADER]=str2num(words{1});
  %
 otherwise
  error(['Unknown entry : ',entry])
end

return