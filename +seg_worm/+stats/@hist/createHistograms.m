function all_hists = createHistograms(feature_file_paths)
%worm2histogram   Convert worm features to their histogram.
%
%   seg_worm.stats.hist.createHistograms
%
%   OLD_NAME: 
%   - seg_worm.w.stats.worm2histogram 
%
%   Inputs
%   =======================================================================
%   feature_file_paths : (char or cellstr), list of feature files to
%   convert into histograms
%
%   Ouputs
%   =======================================================================
%   all_hists : 
%
%
%   TODO: refactor this to be part of the constructor instead of a separate
%   method
%
%   TODO: Could create code that would allow saving in some way without
%   merging or code that would allow returning as a cell array without
%   merging
%
%
%   JAH: Currently working on this file ...
%
%   Code snippet for running in:
%   seg_worm.stats.testing_code
%
%   See Also:
%   seg_worm.stats.hist.initObject

%Loop over all feature files and get histogram objects for each
%--------------------------------------------------------------------------
if ischar(feature_file_paths)
    feature_file_paths = {feature_file_paths};
end

n_files = length(feature_file_paths);
hist_cell_array = cell(1,n_files);

for iFile  = 1:n_files
    hist_cell_array{iFile} = h__getHistsFromFeaturePath(feature_file_paths{iFile});
end

%Merge the objects from each file
%--------------------------------------------------------------------------
all_hists = seg_worm.stats.hist.mergeObjects(hist_cell_array);

end

%==========================================================================
function hist_objs = h__getHistsFromFeaturePath(feature_file_path)

% Initialize the worm data information.
%--------------------------------------------------------------------------
h = load(feature_file_path);

%Movement histograms - DONE
m_hists = h_computeMHists(h,seg_worm.stats.movement_specs.getSpecs);

%Simple histograms - DONE
s_hists = h_computeSHists(h,seg_worm.stats.simple_specs.getSpecs);

%Event histograms - DONE
e_hists = h_computeEHists(h,seg_worm.stats.event_specs.getSpecs);

hist_objs = [m_hists s_hists e_hists];

end

%Converting data to histograms
%==========================================================================
function e_hists = h_computeEHists(h,specs)

n_specs    = length(specs);
temp_hists = cell(1,n_specs);

for iSpec = 1:n_specs
   
   cur_specs = specs(iSpec);

   %NOTE: Because we are doing structure array indexing, we need to capture
   %multiple outputs using [], otherwise we will only get the first value
   %...
   cur_data  = eval(['h.' cur_specs.feature_field]);
    
   if ~isempty(cur_data) && ~isempty(cur_specs.sub_field)
      %This will go from:
      %   frames (structure array) 
      %to:
      %   frames.time
      %for example.
      %
      %It is also used for event.ratio.time and event.ratio.distance
      %going from:
      %ratio (structure or [])
      %to:
      %   ratio.time
      %   ratio.distance
      %
      %
      cur_data = [cur_data.(cur_specs.sub_field)];
   end
   
   temp_hists{iSpec} = seg_worm.stats.hist(cur_data,cur_specs,'event','all','all');
   
end

e_hists = [temp_hists{:}];


end

function m_hists = h_computeMHists(h,specs)
%
%   m_hists = h_computeMHists(h,specs)
%
%   Inputs
%   =======================================================================
%   
%

MAX_NUM_HIST_OBJECTS = 1000;

motion_modes = h.worm.locomotion.motion.mode;

n_frames = length(motion_modes);

indices_by_motion = {1:n_frames ...
    find(motion_modes == 1) ...
    find(motion_modes == 0) ...
    find(motion_modes == -1)};

motion_types = {'all' 'forward' 'paused' 'backward'};
data_types   = {'all' 'absolute' 'positive' 'negative'};

all_hist_objects = cell(1,MAX_NUM_HIST_OBJECTS);
hist_count = 0; 

n_specs = length(specs);

for iSpec = 1:n_specs
   
   cur_specs = specs(iSpec);
   cur_data  = eval(['h.' cur_specs.feature_field]);
    
   if ~isnan(cur_specs.index)
      %This is basically for eigenprojections
      %JAH: I really don't like the orientation: [Dim x n_frames]
      cur_data = cur_data(cur_specs.index,:);
   end
   
   for iMotion = 1:4 
      cur_motion_type = motion_types{iMotion};
       
      hist_count = hist_count + 1;
      temp_data  = cur_data(indices_by_motion{iMotion});
      
      %temp_data(isinf(temp_data) | isnan(temp_data)) = [];
      
      %TODO: This can be sped up significantly as the histogram
      %calculations are redundant ...
      
      all_hist_objects{hist_count} = seg_worm.stats.hist(temp_data,cur_specs,'motion',cur_motion_type,data_types{1});
      if cur_specs.is_signed
         all_hist_objects{hist_count+1} = seg_worm.stats.hist(abs(temp_data),          cur_specs,'motion',cur_motion_type,data_types{2});
         all_hist_objects{hist_count+2} = seg_worm.stats.hist(temp_data(temp_data > 0),cur_specs,'motion',cur_motion_type,data_types{3});
         all_hist_objects{hist_count+3} = seg_worm.stats.hist(temp_data(temp_data < 0),cur_specs,'motion',cur_motion_type,data_types{4});
         hist_count = hist_count + 3;     
      end
   end
end

m_hists = [all_hist_objects{1:hist_count}];

end

function s_hists = h_computeSHists(h,specs)

n_specs = length(specs);

temp_hists = cell(1,n_specs);

for iSpec = 1:n_specs
   
   cur_specs = specs(iSpec);
   cur_data  = eval(['h.' cur_specs.feature_field]);
    
   temp_hists{iSpec} = seg_worm.stats.hist(cur_data,cur_specs,'simple','all','all');
   
end

s_hists = [temp_hists{:}];


end


