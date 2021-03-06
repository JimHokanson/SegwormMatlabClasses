function hist_objs = initObjects(obj,feature_obj)
%
%   hist_objs = seg_worm.stats.hist.manager.initObjects(feature_file_paths)
%
%   INPUTS
%   =======================================
%   feature_obj : (seg_worm.features or strut) This may truly be a feature
%       object or the old structure. Both have the same format.
%
%   This is essentially the constructor code. I moved it in here to avoid
%   the indenting.
%
%   Improvements
%   -----------------------------------------------------------------------
%   - We could optimize the histogram calculation data for motion data
%   - for event data, we should remove partial events (events that start at
%   the first frame or end at the last frame)

%Movement histograms - DONE
m_hists = h_computeMHists(feature_obj,seg_worm.stats.movement_specs.getSpecs);

%Simple histograms - DONE
s_hists = h_computeSHists(feature_obj,seg_worm.stats.simple_specs.getSpecs);

%Event histograms - DONE

%:/ HACK
n_samples = length(feature_obj.morphology.length);

e_hists = h_computeEHists(feature_obj,seg_worm.stats.event_specs.getSpecs,n_samples);

hist_objs = [m_hists s_hists e_hists]';

end

function obj = h__createIndividualObject(data,specs,hist_type,motion_type,data_type)
%
%
%   This method actually creates and populates the object
%
%   Inputs
%   -----------------------------------------------------------------------
%
%
%   TODO:
%   - consider preinstantiating objects
%

obj  = h__getObject(length(data),specs,hist_type,motion_type,data_type);

if obj.n_samples == 0
    return
end

% Compute the histogram counts for all the data.
%-------------------------------------------------------------
[obj.bins,edges] = h__computeBinInfo(data,obj.resolution);

counts = histc(data, edges);
counts(end) = []; %Remove the extra bin at the end (for overflow)

obj.counts = counts;
obj.pdf    = counts./sum(counts);

%Compute stats
%--------------------------------------------------------------
h__computeStats(obj,data)

end

function obj = h__getObject(n_samples,specs,hist_type,motion_type,data_type)
%
%   
%
%

obj = seg_worm.stats.hist;

h__initMeta(obj,specs,hist_type,motion_type,data_type)

obj.n_samples = n_samples;

end

function data = h__filterData(data)
data(isinf(data) | isnan(data)) = [];
end

function mask = h__getFilterMask(data)
%
%   I wasn't sure which of these functions I wanted

mask = isinf(data) | isnan(data);

end

function h__computeStats(obj,data)


obj.mean_per_video  = mean(data);

%I couldn't resist optimizing this since we've already calculated the mean
n_samples = length(data);
if n_samples == 1
    obj.std_per_video = 0;
else
    obj.std_per_video = sqrt(1/(n_samples-1)*sum((data - obj.mean_per_video).^2));
end

end

function h__initMeta(obj,specs,hist_type,motion_type,data_type)
obj.field            = specs.getLongField();
obj.feature_category = specs.feature_category;
obj.resolution       = specs.resolution;
obj.is_zero_bin      = specs.is_zero_bin;
obj.is_signed        = specs.is_signed;
obj.name             = specs.name;
obj.short_name       = specs.short_name;
obj.units            = specs.units;

obj.hist_type   = hist_type;
obj.motion_type = motion_type;
obj.data_type   = data_type;
end

function [bins,edges] = h__computeBinInfo(data,resolution)
%
%
%   NOTE: This version may have an extra bin than the previous version but
%   this one is MUCH simpler and merging should be much simpler as edges
%   should always align ...
%
%
%   min -65.4
%   max 20.01
%   resolution 1
%   Old:
%   edges -65.5 to 20.5
%   New:
%   edges -70 to 21
%
%


MAX_NUMBER_BINS = 1e6; %The maximum # of bins that we'll use. Since the data
%is somewhat random, outliers could really chew up memory. I'd prefer not
%to have some event which all of a sudden tells the computer we need to
%allocate a few hundred gigabytes of data. If this does ever end up a
%problem we'll need a better solution (or really A solution)

%Compute the data range & padding
%------------------------------------------------
min_data = min(data);
max_data = max(data);

min_edge = floor(min_data/resolution)*resolution;
max_edge = ceil(max_data/resolution)*resolution;

%If we have a singular value, then we will get a singular edge, which isn't
%good for binning. We always need to make sure that our data is bounded by
%a high and low end. Given how hist works (it is inclusive on the low end,
%when we only have one edge we add a second edge by increasing the high
%end, NOT by decreasing the low end.
%
%i.e. in Matlab you can't bound 3 by having edges at 2 & 3, the edges would
%need to be at 3 & 4

if min_edge == max_edge
    max_edge = min_edge + resolution;
end

n_values = (max_edge - min_edge)/resolution + 1;

if n_values > MAX_NUMBER_BINS
    %TODO: Make the error more explicit
    error('Given specified resolution there are too many data bins')
end

edges = min_edge:resolution:max_edge;
bins  = edges(1:end-1) + resolution/2;

end

%Converting data to histograms
%==========================================================================
function e_hists = h_computeEHists(h,specs,n_samples)

n_specs    = length(specs);
temp_hists = cell(1,4*n_specs);

cur_hist_index = 0;

for iSpec = 1:n_specs
    cur_specs = specs(iSpec);
    
    %seg_worm.stats.event_specs.getData
    cur_data = cur_specs.getData(h,n_samples);
    
    cur_data = h__filterData(cur_data);
    
    if cur_specs.is_signed
        temp_hists{cur_hist_index+1} = h__createIndividualObject(cur_data,cur_specs,'event','all','all');
        temp_hists{cur_hist_index+2} = h__createIndividualObject(abs(cur_data),cur_specs,'event','all','absolute');
        positive_mask = cur_data > 0;
        negative_mask = cur_data < 0;
        temp_hists{cur_hist_index+3} = h__createIndividualObject(cur_data(positive_mask),cur_specs,'event','all','positive');
        temp_hists{cur_hist_index+4} = h__createIndividualObject(cur_data(negative_mask),cur_specs,'event','all','negative');
        cur_hist_index    = cur_hist_index + 4;
    else
        cur_hist_index    = cur_hist_index + 1;
        temp_hists{cur_hist_index} = h__createIndividualObject(cur_data,cur_specs,'event','all','all');
    end
    
end

e_hists = [temp_hists{1:cur_hist_index}];


end

function m_hists = h_computeMHists(feature_obj,specs)
%
%   m_hists = h_computeMHists(h,specs)
%
%   For movement features, we compute either 4 or 16 histogram objects,
%   depending on whether or not the feature can be signed. If the data is
%   not signed, then we compute 4. If it is, we compute 16. We compute
%   histograms for when the motion of the midbody is:
%       - going forward
%       - going backward
%       - paused
%       - all 3 of the above combined
%
%   If signed, this also gets computed on feature data that is:
%       - positive
%       - negative
%       - absolute value of data
%       - both positive and negative
%
%   By combining the motion of the midbody with the sign of the data, we
%   get 16 different possible combinations
%
%
%   Inputs
%   =======================================================================
%   h :
%   specs : ( )
%
%   Improvements:
%   =======================================================================
%   - we could significantly reduce the amount of binning done in this
%   function

MAX_NUM_HIST_OBJECTS = 1000;

%---------------------------------------------------------
motion_modes = feature_obj.locomotion.motion.mode;

n_frames = length(motion_modes);

indices_use_mask = {...
    true(1,n_frames) ...
    motion_modes == 1 ...
    motion_modes == 0 ...
    motion_modes == -1};

%NOTE: motion types refers to the motion of the worm's midbody
motion_types = {'all' 'forward'     'paused'    'backward'};
data_types   = {'all' 'absolute'    'positive'  'negative'};
%---------------------------------------------------------

all_hist_objects = cell(1,MAX_NUM_HIST_OBJECTS);
hist_count = 0;

n_specs = length(specs);

for iSpec = 1:n_specs
    
    cur_specs = specs(iSpec);
    
    cur_data = cur_specs.getData(feature_obj);
    
    good_data_mask = ~h__getFilterMask(cur_data);
    
    for iMotion = 1:4
        cur_motion_type = motion_types{iMotion};
        
        hist_count = hist_count + 1;
        temp_data  = cur_data(indices_use_mask{iMotion} & good_data_mask);
        
        all_obj = h__createIndividualObject(temp_data,cur_specs,'motion',cur_motion_type,data_types{1});
        all_hist_objects{hist_count} = all_obj;
        if cur_specs.is_signed
            
            %TODO: This could be improved by merging results from positive and
            %negative ...
            all_hist_objects{hist_count+1} = h__createIndividualObject(abs(temp_data),cur_specs,'motion',cur_motion_type,data_types{2});
            
            
            %NOTE: To get a speed up, we don't rely on h__createIndividualObject
            %Instead we take the positive and negative aspects of the
            %object that included all data.
            
            %positive object ----------------------------------------
            pos_obj  = h__getObject(0,cur_specs,'motion',cur_motion_type,data_types{3});
            
            I_pos = find(all_obj.bins > 0 & all_obj.counts > 0,1);
            
            if ~isempty(I_pos)
                pos_obj.bins      = all_obj.bins(I_pos:end);
                pos_obj.counts    = all_obj.counts(I_pos:end);
                pos_obj.n_samples = sum(pos_obj.counts);
                
                h__computeStats(pos_obj,temp_data(temp_data > 0))
            end
            
            %negative object -----------------------------------------
            neg_obj  = h__getObject(0,cur_specs,'motion',cur_motion_type,data_types{4});
            
            I_neg = find(all_obj.bins < 0 & all_obj.counts > 0,1,'last');
            
            if ~isempty(I_neg)
                neg_obj.bins      = all_obj.bins(1:I_neg);
                neg_obj.counts    = all_obj.counts(1:I_neg);
                neg_obj.n_samples = sum(neg_obj.counts);
                h__computeStats(neg_obj,temp_data(temp_data < 0))
            end
            
            
            
            %final assignments --------------------------------------
            all_hist_objects{hist_count+2} = pos_obj;
            all_hist_objects{hist_count+3} = neg_obj;
            hist_count = hist_count + 3;
        end
    end
end

m_hists = [all_hist_objects{1:hist_count}];

end

function s_hists = h_computeSHists(feature_obj,specs)

n_specs = length(specs);

temp_hists = cell(1,n_specs);

for iSpec = 1:n_specs
    cur_specs = specs(iSpec);
    cur_data  = h__filterData(cur_specs.getData(feature_obj));
    
    temp_hists{iSpec} = h__createIndividualObject(cur_data,cur_specs,'simple','all','all');
end

s_hists = [temp_hists{:}];


end
