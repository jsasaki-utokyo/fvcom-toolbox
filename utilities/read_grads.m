function [data,header]=read_grads(file_name,var_name,varargin)
% [var,header]=read_grads(file_name,var_name) is a function to read
% in binary files using their GrADS descriptor header. 
% If var_name is 'all' or non-existent, all variables are read and 
% sent to the base-workspace.
% Please check before-hand the data-type (big/little endian), as well
% as the data precision (single/double...), and adjust if necessary 
% the header at the beginning of read_grads.m.
% Optionally, you can improve the indications from the header file in 
% the sub-function grads_name (hereafter).
% 
% Kristof Sturm, sturm@dkrz.de or sturm@lgge.obs.ujf-grenoble.fr, 29.10.02

% Overview
% I introduce here an interface between GrADS and Matlab: the read_grads.m 
% routine reads the GrADS control file and retrieves the variable from the 
% binary file accordingly. As far as possible, the components of the GrADS 
% Data Descriptor File (or control file) were implemented in this routine. 
% The companion routine write_grads.m writes a variable in the Matlab 
% workspace as a binary file, with the appropriate control file. 
% Examples are given hereafter.
% Example:
% [data,header]=read_grads('filename.ctl',''); % Reads only the header
% [data,header]=read_grads('filename.ctl','varname'); % Reads the header and <varname>
% [data,header]=read_grads('filename.ctl','varname','z',[z1,z2],'lon',[lon1,lon2],'lat',[lat1,lat2],'t',[t1,t2]);
% [data,header]=read_grads('filename.nc','varname'); % Extension of read_grads.m to NetCDF variables
% * data is a 4-D variable, with X/lon as the first dimension, 
% Y/lat as second dimension, Z/lev as third dimension and time in fourth dimension.
% * header is a structure variable, which contains all relevant information 
% from the GrADS control file. In particular, the header.XDEF and header.YDEF 
% sub-structures contain information for plotting the variable. header.XDEF.num 
% is the number of X fields, and header.XDEF.vec the X coordinates.
% Further details can be seen from the dedicated page: http://www.uib.no/People/kst064/REMOiso-toolbox/Matlab_Toolbox.html

global l_quiet

n=strmatch('quiet',{varargin{1:2:end}});
if isempty(n)
  l_quiet=0;
else
  switch lower(varargin{2*(n-1)+2})
   case {1,'yes'}
    l_quiet=1;
   case {0,'no'}
    l_quiet=0;
  end
end

header=struct('FILENAME',file_name,...
	      'VARSIZE',{[]},...
	      'NVAR',0,...
	      'DATANAME','-',...
	      'FID',{[]},...
	      'BINTYPE','-',...
	      'BINPRECISION','float32',...
	      'DSET','-',...
	      'DTYPE',1,...
	      'INDEX','-',...
	      'TITLE','-',...
	      'UNDEF',0,...
	      'OPTIONS',{[]},...
	      'XDEF',{[]},...
	      'YDEF',{[]},...
	      'ZDEF',{[]},...
	      'TDEF',{[]},...
	      'VARS',struct([]),...
	      'FILEHEADER',0,...
	      'THEADER',0,...
	      'XYHEADER',0);

header.XDEF.rev=0;
header.YDEF.rev=0;
header.ZDEF.rev=0;

% [comp,max_size,endian]=computer;
endian='L';
if endian=='L'
  header.BINTYPE='ieee-le'; 
elseif endian=='B'
  header.BINTYPE='ieee-be';
end

data=[];

file_name=deblank(file_name);
if exist(file_name,'file')
  fid_ctl = fopen(file_name,'r');
  % disp(['CTL opened: ',num2str(fid_ctl)])
  % [header.FID]=deal([fid_ctl 0]);
  header(1).FID{1}=fid_ctl;
  if isunix, sep='/';
  else
    sep='\\';
  end
  if ~isempty(strfind(file_name,sep))
    slash=max(strfind(file_name,sep));
    header(1).DIR=file_name(1:slash);
    file_name=file_name((slash+1):end);
  end
  header(1).FILENAME=file_name;
else
  error(['Inexistent header file ',file_name])
end

% Enabling read_grads to handle netcdf files. 

point=strfind(file_name,'.');
suf=lower(file_name(point(end)+1:end));

if strcmp(suf,'ctl')
% reading the .ctl file : 
    while ~feof(fid_ctl)
      line=fgets(fid_ctl);
      if isempty(line(1:end-1))
          break
      end
      % cutting the line into words.
      line_gaps=isspace(line);
      % removing consecutive blancks.
      while line_gaps(1)==1
          line_gaps(1)=[];
          line(1)=[]; 
      end
      if line(1)=='*'
        if ~l_quiet
          disp(line)
        end
      else
        no_gaps=0;
        while ~no_gaps
          line_gap2=find(line_gaps);
          if any(diff(line_gap2)==1)
        dup=1+find(diff(line_gap2)==1);
        line_gaps(line_gap2(dup))=[];
        line(line_gap2(dup))=[];
          else
        no_gaps=1;
          end
        end
        entry=line(1:(min(find(line_gaps))-1));
        words={''};
        for n=1:(sum(line_gaps)-1)
          words{n}=line((line_gap2(n)+1):(line_gap2(n+1)-1));
        end
        header=read_header(entry,words);
        % disp(words)
      end
    end

% Opening the binary file for later reading:
fid_bin = fopen(header.DATANAME,'r', header.BINTYPE );
% disp(['BIN opened: ',num2str(fid_bin)])
[header(1).FID] = {fid_ctl, fid_bin};

fclose(fid_ctl);
nvar=size(header.VARS,1);
header.NVAR=nvar;

if nargin>1 & isempty(var_name)
  if ~l_quiet
    disp(['Header for ',header.FILENAME,' read.'])
  end
  fclose(header.FID{2});
  return
end

if nargout==0
  assignin('caller','header',header)
  assignin('caller','vars',vars);
end

if header.DTYPE          % Gridded data set
  
  % Computing the size of each variable:
  var_size=zeros(nvar,4);  % array of the var. size [X,Y,Z,T]
  var_size(:,1)=header.XDEF.num;
  var_size(:,2)=header.YDEF.num;
  var_size(:,3)=cat(1,header.VARS.levs);
  var_size(:,4)=header.TDEF.num;
  header.VARSIZE={var_size};
  if exist('grads_name','file')
    [vars,header]=grads_name(header,vars);
  end
  
  % reading the bin. file :
  if nargin == 1 | strcmp(lower(var_name),'all')
    read_data(header);
  elseif nargin == 2
    t_limits=[1 header.TDEF.num];
    ivar=strmatch(var_name,{header.VARS.name},'exact');
    if isempty(ivar)
      if ~l_quiet
	disp(['The variable ',var_name,' cannot be found in ',header.DATANAME])
      end
      data=[];
      return
    end
    z_limits=[1 header.VARS(ivar).levs];
    y_limits=[1 header.YDEF.num];
    x_limits=[1 header.XDEF.num];
    data=read_var(header,var_name,t_limits,z_limits,y_limits,x_limits);
    header.XDEF.vec=header.XDEF.vec(x_limits(1):x_limits(2));
    header.XDEF.num=length(header.XDEF.vec);
    header.YDEF.vec=header.YDEF.vec(y_limits(1):y_limits(2));
    header.YDEF.num=length(header.YDEF.vec);
  else
    for n=1:2:nargin-2
      switch lower(varargin{n})
       case {'x'}
	x_limits=varargin{n+1};
       case {'lon'}
	lon_limits=sort(varargin{n+1});
	x_limits=interp1(header.XDEF.vec,1:header.XDEF.num,...
			 lon_limits,'nearest',NaN);
	if any(isnan(x_limits))
	  if ~l_quiet
	    disp(['The X-limits ',num2str(lon_limits),...
		  ' exceed the data coverage ',...
		  num2str(header.XDEF.vec([1 end]))])
	  end
	  data=[];
	  return
	end
       case {'y'}
	y_limits=varargin{n+1};
       case {'lat'}
	lat_limits=sort(varargin{n+1});
	y_limits=interp1(header.YDEF.vec,1:header.YDEF.num,...
			 lat_limits,'nearest',NaN);
	if any(isnan(y_limits))
	  if ~l_quiet
	    disp(['The Y-limits ',num2str(lat_limits),...
		  ' exceed the data coverage ',...
		  num2str(header.YDEF.vec([1 end]))])
	  end
	  data=[];
	  return
	end
	if header.YDEF.rev,
	  y_limits=header.YDEF.num*ones(size(y_limits))-fliplr(y_limits)+1;
	end
       case {'t','time'}
	t_limits=varargin{n+1};
       case {'z','lev'}
	z_limits=varargin{n+1};
	ivar=strmatch(var_name,{header.VARS.name},'exact');
	if isempty(ivar)
	  if ~l_quiet
	    disp(['The variable ',var_name,' cannot be found in ',header.DATANAME])
	  end
	  data=[];
	  return
	end
	if any(z_limits>header.VARS(ivar).levs)
	  if ~l_quiet
	    disp(['The variable ',var_name,' has only ',...
		  header.VARS(ivar).levs,' levels.'])
	  end
	  data=[];
	  return
	end
      end
    end
    if ~exist('t_limits','var'), t_limits=[1 header.TDEF.num]; end
    if ~exist('z_limits','var')
      ivar=strmatch(var_name,{header.VARS.name},'exact');
      if isempty(ivar)
	if ~l_quiet
	  disp(['The variable ',var_name,' cannot be found in ',header.DATANAME])
	end
	data=[];
	return
      end
      z_limits=[1 header.VARS(ivar).levs];
    end
    if ~exist('y_limits','var')
      y_limits=[1 header.YDEF.num];
    end
    if ~exist('x_limits','var')
      x_limits=[1 header.XDEF.num];
    end      
    data=read_var(header,var_name,t_limits,z_limits,y_limits,x_limits);
    header.XDEF.vec=header.XDEF.vec(x_limits(1):x_limits(2));
    header.XDEF.num=length(header.XDEF.vec);
    if header.YDEF.rev,
      y_limits=header.YDEF.num*ones(size(y_limits))-fliplr(y_limits)+1;
    end
    header.YDEF.vec=header.YDEF.vec(y_limits(1):y_limits(2));
    header.YDEF.num=length(header.YDEF.vec);
  end
  
  %    if header.XDEF.rev, header.XDEF.vec=header.XDEF.vec([end:-1:1]); end
  %    if header.YDEF.rev, header.YDEF.vec=header.YDEF.vec([end:-1:1]); end
  %    if header.ZDEF.rev, header.ZDEF.vec=header.ZDEF.vec([end:-1:1]); end
  
else                     % Station data set
  
  data=read_station(header);
  
end

fclose(header.FID{2});

elseif strcmp(suf,'nc') % Enabling NetCDF in read_grads, based on snctools
  old_dir=pwd;
  if isfield(header,'DIR')
    cd(header.DIR)
  end % if
  
  % Reading the header, conforming to GrADS conventions:
  ncheader=nc_info(file_name);
  
  header.FILENAME=ncheader.Filename;
  header.DATANAME=ncheader.Filename;
  header.NVAR=length(ncheader.DataSet);

  fields={{'lon','longitude','x','i'},'XDEF';
          {'lat','latitude','y','j'},'YDEF';
          {'lev','plev','nv','bnds'},'ZDEF';
          {'time','tps','t','l'},'TDEF'};
  for idim=1:size(fields,1)
    for dim=fields{idim,1}
      dim_id=strmatch(dim{1},{ncheader.DataSet.Name},'exact');
      if ~isempty(dim_id), break, end
    end % for
    if isempty(dim_id)
      disp(['Dim ',fields{idim,2},' is missing.'])
      vec=1;
      num=1;
      type='LINEAR';
    else
      num=ncheader.DataSet(dim_id).Size;
      vec=nc_varget(file_name,dim{1});
      if length(unique(diff(vec)))==1
        type='LINEAR';
      else
        type='LEVELS';
      end
      if all(diff(vec)<0)
        rev=1;
        vec=flipud(reshape(vec,[],1));
      elseif all(diff(vec)>0)
        rev=0;
      else
        error(['Unsuitable vector ',dim{1}])
      end
    end % if
    
    DEF=struct('num',num,'vec',vec,'type',type,'rev',rev);
    
    header=setfield(header,fields{idim,2},DEF);
  end % for idim

  for ivar=1:header.NVAR
    header.VARS(ivar).id=ivar;
    header.VARS(ivar).name=ncheader.DataSet(ivar).Name;
    header.VARS(ivar).size=ones(1,4);
    for idim=1:length(ncheader.DataSet(ivar).Dimension)
      [dim_id,j]=ind2sub(size(cat(1,fields{:,1})),...
                         strmatch(ncheader.DataSet(ivar).Dimension(idim),...
                             cat(1,fields{:,1}),'exact'));
      header.VARS(ivar).size(dim_id)=...
          getfield(eval(['header.',fields{dim_id,2}]),'num');
    end % for
    header.VARS(ivar).levs=header.VARS(ivar).size(3);
  end % for

  header.VARSIZE={cat(1,header.VARS.size)};
  
  if nargin==1 | isempty(var_name) | strcmp(var_name,'all')
    data=[];
    if exist('old_dir','var')
      cd(old_dir)
    end
    return
  end
  
  var_id=strmatch(var_name,{ncheader.DataSet.Name});
  if isempty(var_id)
    error(['Variable ',var_name,' is not available in ',file_name,...
           ': ',{ncheader.DataSet.Name}])
  end % if 
  
  for n=1:2:nargin-2
    switch lower(varargin{n})
     case {'x'}
      x_limits=varargin{n+1};
     case {'lon'}
      lon_limits=sort(varargin{n+1});
      x_limits=interp1(header.XDEF.vec,1:header.XDEF.num,...
                       lon_limits,'nearest',NaN);
      if any(isnan(x_limits))
        if ~l_quiet
          disp(['The X-limits ',num2str(lon_limits),...
                ' exceed the data coverage ',...
                num2str(header.XDEF.vec([1 end])')])
        end
        data=[];
        return
      end
     case {'y'}
      y_limits=varargin{n+1};
     case {'lat'}
      lat_limits=sort(varargin{n+1});
      y_limits=interp1(header.YDEF.vec,1:header.YDEF.num,...
                       lat_limits,'nearest',NaN);
      if any(isnan(y_limits))
        if ~l_quiet
          disp(['The Y-limits ',num2str(lat_limits),...
                ' exceed the data coverage ',...
                num2str(header.YDEF.vec([1 end])')])
        end
        data=[];
        return
      end
      if header.YDEF.rev,
        y_limits=header.YDEF.num*ones(size(y_limits))-fliplr(y_limits)+1;
      end
     case {'t','l','time'}
      t_limits=varargin{n+1};
     case {'z','k','lev'}
      z_limits=varargin{n+1};
      ivar=strmatch(var_name,{header.VARS.name},'exact');
      if isempty(ivar)
        if ~l_quiet
          disp(['The variable ',var_name,' cannot be found in ',header.DATANAME])
        end
        data=[];
        return
      end
      if any(z_limits>header.VARS(ivar).levs)
        if ~l_quiet
          disp(['The variable ',var_name,' has only ',...
                header.VARS(ivar).levs,' levels.'])
        end
        data=[];
        return
      end
    end % switch
  end % for
  
  if ~exist('t_limits','var'), t_limits=[1 header.TDEF.num]; end
  if ~exist('z_limits','var')
    ivar=strmatch(var_name,{header.VARS.name},'exact');
    if isempty(ivar)
      if ~l_quiet
        disp(['The variable ',var_name,' cannot be found in ',header.DATANAME])
      end
      data=[];
      return
    end
    z_limits=[1 header.VARS(ivar).levs];
  end % if
  if ~exist('y_limits','var')
    y_limits=[1 header.YDEF.num];
  end
  if ~exist('x_limits','var')
    x_limits=[1 header.XDEF.num];
  end      
  % data=read_var(header,var_name,t_limits,z_limits,y_limits,x_limits);
  % Retrieve the corresponding variable:
  
  if all(header.VARS(var_id).size([3 4])==[1 1])
    data=nc_varget(file_name,var_name,...
                   [y_limits(1),x_limits(1)]-1,...
                   diff(cat(1,y_limits,x_limits),1,2)'+1);
  elseif header.VARS(var_id).size(3)==1
    data=nc_varget(file_name,var_name,...
                   [t_limits(1),y_limits(1),x_limits(1)]-1,...
                   diff(cat(1,t_limits,y_limits,x_limits),1,2)'+1);
  elseif header.VARS(var_id).size(4)==1
    data=nc_varget(file_name,var_name,...
                   [z_limits(1),y_limits(1),x_limits(1)]-1,...
                   diff(cat(1,z_limits,y_limits,x_limits),1,2)'+1);
  else
    data=nc_varget(file_name,var_name,...
                   [t_limits(1),z_limits(1),y_limits(1),x_limits(1)]-1,...
                   diff(cat(1,t_limits,z_limits,y_limits,x_limits),1,2)'+1);
  end
  
  header.XDEF.vec=header.XDEF.vec(x_limits(1):x_limits(2));
  header.XDEF.num=length(header.XDEF.vec);
  if header.YDEF.rev,
    y_limits=header.YDEF.num*ones(size(y_limits))-fliplr(y_limits)+1;
  end 
  if diff(y_limits)<0, error('Unsuitable Y limits'), end
  header.YDEF.vec=header.YDEF.vec(y_limits(1):y_limits(2));
  header.YDEF.num=length(header.YDEF.vec);
  
  header.ZDEF.vec=z_limits(1):z_limits(2);
  header.ZDEF.num=diff(z_limits)+1;
  
  header.TDEF.vec=t_limits(1):t_limits(2);
  header.TDEF.num=diff(t_limits)+1;
  
  if length(size(data))==2
    data=permute(data,[2 1]);
  elseif length(size(data))==3 & header.ZDEF.num==1
    data=permute(data,[3 2 4 1]);
  elseif length(size(data))==3 & header.TDEF.num==1
    data=permute(data,[3 2 1]);
  elseif length(size(data))==4
    data=permute(data,[4 3 2 1]);
  end % if
  
  if header.YDEF.rev
    data=flipdim(data,2);
  end
  
  if exist('old_dir','var')
    cd(old_dir)
  end
else
  disp(['Suffix .',suf,' is not yet enabled.'])
end % if
return
