tt=sol{6};
rks=tt_ranks(tt);
rks=[1,rks',1];
%Go through all the tensor and try to merge modes, if they are amenable to
%merge

tt0=tt;
s=1;
possible_to_merge=true;
while ( possible_to_merge )
  rks = tt_ranks(tt0);
  rks=[1,rks',1];
  sz=tt_size(tt0);
  d=numel(sz);
  %Find first possible dimension that can be merged
  possible_to_merge=false;
  s=1;
  while ( ~possible_to_merge && s < d )
     %Check certain (simple) inequality
     %merging modes s & s+1
     old_mem=sz(s)*rks(s)*rks(s+1)+sz(s+1)*rks(s+1)*rks(s+2);
     new_mem=rks(s)*sz(s)*sz(s+1)*rks(s+2);
     if ( new_mem < old_mem ) 
        tt0=tt_merge(tt0,s);
        possible_to_merge=true;
     end
     s=s+1;
  end
end