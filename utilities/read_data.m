function data=read_data(header,var_name)
% data=read_data(header,var_name) is a function, usually called by
% read_grads.m to use the information contained in the header file
% (.ctl) for reading a binary GrADS file.
% If var_name is 'all', all variables are read and sent to the base-
% workspace.  
% Kristof Sturm, 29.10.02

global l_quiet

vars=header.VARS;
nvar=header.NVAR;
var_size=header.VARSIZE{:};

if nargin == 1 | strcmp(lower(var_name),'all')           % All variables loaded into workspace
  data_list={};
  % initialising the data matrices : 
  var_count=zeros(nvar,1);
  for ivar=1:nvar
    eval([vars(ivar).name,'=NaN*ones(',int2str(header.XDEF.num),',',...
	  int2str(header.YDEF.num),',',...
	  int2str(vars(ivar).levs),',',...
	  int2str(header.TDEF.num),');']);
    fprintf('Size of %s : [%s].\n            %s[X   Y   Z    T]\n',...
	    vars(ivar).name,int2str(size(eval(vars(ivar).name))),...
	    char(32*ones(1,length(vars(ivar).name))));
  end
  
  method = 'loop';  % mainly for debugging... 
    
  switch method
   case 'loop'
    frewind(header.FID{2});
    % Implementing the sequential option:
    if isfield(header,'SEQ') && header.SEQ
      fseek(header.FID{2},4,'bof');
    end

    if header.FILEHEADER
      fseek(header.FID{2},header.FILEHEADER,'bof');
% $$$ if isfield(header,'SEQ') && header.SEQ
% $$$ 	fseek(header.FID{2},8,'cof');
% $$$ end
    end
    
    for l=1:header.TDEF.num
      if header.THEADER
	fseek(header.FID{2},header.THEADER,'cof');
% $$$ 	if isfield(header,'SEQ') && header.SEQ
% $$$ 	  fseek(header.FID{2},8,'cof');
% $$$ 	end
      end
      for ivar=1:nvar
	data=eval(vars(ivar).name);
	count0=var_count(ivar);
	for k=1:vars(ivar).levs
	  if header.XYHEADER
	    fseek(header.FID{2},header.XYHEADER,'cof');
	  end
	  [data(:,:,k,l) count1]=fread(header.FID{2},[header.XDEF.num header.YDEF.num],header.BINPRECISION);
	  count0=count0+count1;
	  if isfield(header,'SEQ') && header.SEQ
	    fseek(header.FID{2},8,'cof');
	  end
	end                                   % end k-loop
	data(data==header.UNDEF)=NaN;   % replacing the
                                        % missing/undefined values
                                        % by NaN (not-a-number)
					
%	data(data==header.UNDEF)=0;   % replacing the
                                      % missing/undefined values by
                                      % NaN (not-a-number)
%	data=sparse(data);			      
	eval([vars(ivar).name,'=data;']);
	var_count(ivar)=var_count(ivar)+count0;
	% bytes(ivar)=ftell(header.FID{2});
	clear data
      end                                     % end ivar-loop
    end                                       % end l-loop
    
    % Send the relevent variables into the base workspace:
    %	assignin('base','var_count',var_count);
    %       assignin('base','bytes',bytes);
    for ivar=1:nvar
      data=eval(vars(ivar).name);
      if header.XDEF.rev, data=flipdim(data,1); end
      if header.YDEF.rev, data=flipdim(data,2); end
      if header.ZDEF.rev, data=flipdim(data,3); end
      assignin('base',vars(ivar).name,data);
    end
    %                                             % end 'loop'
   otherwise
    for ivar=1:nvar
      data=read_var(header,vars(ivar).name);
      assignin('base',vars(ivar).name,data);
      data_list=cat(1,data_list,{data});
    end
    data=[];
    assignin('caller','data_list',data_list);
  end
else                                          % Individual var. reading
  switch var_name
   case {vars.name}
    fprintf('Variable %s extracted from %s. \n',var_name, ...
	    header.DATANAME);
    data=read_var(header,var_name);
   otherwise
    error(['Variable ',var_name,' not available in ', ...
	   header.DATANAME]);
  end
end

return