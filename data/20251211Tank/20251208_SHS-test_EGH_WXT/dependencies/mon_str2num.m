
function [mnums ]  = mon_str2num(mstrs)
  % return month numbers from month names where mstrs is an array of month names.
  % n galbraith 2010/10/27
  mymstrs = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
  for ii=1:12,
      mymsum=sum(mymstrs{ii},2);
      mysumarr(mymsum) =ii;
  end
  thesums=sum(mstrs,2);
  mnums=ones(length(mstrs),1) * NaN;
  goodms=find(thesums <= length(mysumarr));
  mnums(goodms) = mysumarr(thesums(goodms))';
end
