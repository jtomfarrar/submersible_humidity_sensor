function [newstruct] = struct_trunc(oldstruct, indx)
%
% truncate (or subset) a structure
% usage [newstruct] = struct_trunc(oldstruct, indx)
% where oldstruct is the existing struct, indx is the index to good values
% newstruct contains all the fields in oldstruct, but those fields with
% one dimension that matches the length of mday will be indexed on indx
%
debug=0;
% use mday to figure out where to truncate 2d arrays
ntpts=length(oldstruct.mday);

if debug ==1,
  fprintf('length mday %d length indx %d\n', ...
    ntpts, length(indx));
end

% get length of time param, for 2d arrays later
fno=fieldnames(oldstruct);
for i=1:length(fno)
    field_name=fno{i};
    if debug ==1,
       fprintf('%s ',field_name); end
    cmd=['tempvar =  oldstruct.' field_name ';' ];
    eval (cmd);
    [m,n]=size(tempvar);
    if debug ==1,
      fprintf(' size %d %d\n',n,m);  end
    if m == ntpts,
        % time is first dimension
        cmd=['newstruct.' field_name ' =  tempvar(indx,:);' ];
        eval(cmd);
    else
        if n == ntpts,
            % time is second dimension
            cmd=['newstruct.' field_name ' =  tempvar(:,indx);' ];
            eval(cmd);
        else
            % if no dimension matches time, don't sort
            if debug ==1, fprintf('%s no trunc\n', field_name); end
            cmd=['newstruct.' field_name ' =  tempvar;' ];
            eval(cmd);
        end
    end
end
% removed this: don't sort unless one dim matches mday...
% if m ==1 | n ==1,
% this is a 1-d array
%       cmd=['newstruct.' field_name ' =  tetmpvar(indx);' ];
%       eval(cmd);
