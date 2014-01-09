function initObject(obj,exp_hists,ctl_hists)
%
%
%   seg_worm.stats.manager.initObject(obj,exp_hists,ctl_hists)
%


n_objs = length(exp_hists);

stats_objs(n_objs) = seg_worm.stats();

%:/ Sadly this needs to be done before hand to be the same ...
%It might be updated during object initialization ...
%
%TODO: This looks nicely vectorized, but it breaks the organization
%significantly ...
%
%How much of an impact do we get if we move this to being computed 
%for each object, instead of all of them at once?
p_t_all          = mattest([exp_hists.mean_per_video]',[ctl_hists.mean_per_video]');
[stats_objs.p_t] = sl.struct.dealArray(p_t_all);

for iObj = 1:n_objs
   %seg_worm.stats.initObject
   stats_objs(iObj).initObject(exp_hists(iObj),ctl_hists(iObj)); 
end

obj.stats = stats_objs;

[~,q_t_all]      = mafdr([stats_objs.p_t]);
[stats_objs.q_t] = sl.struct.dealArray(q_t_all);

[~,q_w_all]      = mafdr([stats_objs.p_w]);
[stats_objs.q_w] = sl.struct.dealArray(q_w_all);


obj.p_worm = min([stats_objs.p_w]);
obj.q_worm = min([stats_objs.q_w]);

end