function initObject(obj,exp_hist,ctl_hist)
%worm2StatsInfo  Compute worm statistics information and save it to a file.
%
%   seg_worm.stats.initObject
%
%   See Also:
%   seg_worm.stats.helpers.swtest
%

ALPHA = 0.05; %Not really used since we don't examine H, just p
TAIL  = 0;

obj.field              = exp_hist.field;
obj.name               = exp_hist.name; 
obj.short_name         = exp_hist.short_name;
obj.units              = exp_hist.units;
obj.feature_category   = exp_hist.feature_category;
obj.hist_type          = exp_hist.hist_type;
obj.motion_type        = exp_hist.motion_type;
obj.data_type          = exp_hist.data_type;


is_exclusive = (exp_hist.none_valid && ctl_hist.all_valid) ...
        || (exp_hist.all_valid && ctl_hist.none_valid);

%zscore
%--------------------------------------------------------------------------
%This definition is slightly different than the old version, but matches
%the textual description (in code, what about in published paper?)
if isnan(exp_hist.mean)
    if ctl_hist.n_valid_measurements > 1
        obj.z_score_experiment = -Inf;
    else
        obj.z_score_experiment = NaN;
    end
elseif isnan(ctl_hist.mean)
    if exp_hist.n_valid_measurements > 1
        obj.z_score_experiment = Inf;
    else
        obj.z_score_experiment = NaN;
    end
else
    %This might need to be means_per_video, not the mean ...
    obj.z_score_experiment = (exp_hist.mean - ctl_hist.mean)/ctl_hist.std;
end

%TODO: Move this to the histogram, not here ... These are properties of how
%normal the distributions of the histograms are, and have nothing to do
%with the comparative statistics between the two groups
%------------------------------------------------------------------------
p_fields  = {'p_normal_experiment' 'p_normal_control'};
hist_objs = {exp_hist ctl_hist};

for iObj = 1:2
    cur_field    = p_fields{iObj};
    cur_hist_obj = hist_objs{iObj};
    if cur_hist_obj.n_valid_measurements < 3
        obj.(cur_field) = NaN;
    else
        obj.(cur_field) = seg_worm.stats.helpers.swtest(cur_hist_obj.mean_per_video,ALPHA,TAIL);
    end
end



%Rules are:
%--------------------------------------
%p_t
%
%- not in one, but all in the other - use fexact (Fishers Exact)
%- otherwise use mattest
%
%p_w
%- not in one, but all in the other - use fexact
%- partial in both - use Wilcoxon rank-sum test
%- if in both, set to NaN

if is_exclusive
    %This is a literal trasnlation of the code (I think)
    %I'm a bit confused by it ...
    n_expt  = exp_hist.n_videos;
    n_total = n_expt + ctl_hist.n_videos;
    obj.p_w = seg_worm.stats.helpers.fexact(n_expt,n_total,n_expt,n_expt);
    obj.p_t = obj.p_w;
elseif ~(exp_hist.none_valid || ctl_hist.none_valid)
    %We need a few valid values from both ...
    obj.p_w = ranksum(exp_hist.valid_means,ctl_hist.valid_means);
end

%NOTE: This code is for an invidual object, the corrections are done in the
%manager which is aware of all objects ...


%pWValues - these seem to be the real statistics used ...
%- exclusive - fexact  seg_worm.stats.helpers.fexact
%- ranksum

