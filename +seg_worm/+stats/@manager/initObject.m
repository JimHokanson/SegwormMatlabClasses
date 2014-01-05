function initObject(obj,exp_hists,ctl_hists)

% Find exclusive features.
isNaNWorm    = isnan([wormData.dataMeans]);
isNaNControl = isnan([controlData.dataMeans]);
isExclusive  = (all(isNaNWorm) & ~any(isNaNControl)) | (~any(isNaNWorm) & all(isNaNControl));
exclusiveI   = find(isExclusive);
 
n_objs = length(exp_hists);

stats_objs(n_objs) = seg_worm.stats();

for iObj = 1:n_objs
   stats_objs(iObj).initObject(exp_hists(iObj),ctl_hists(iObj)); 
end

end