function [new_data,good_indices,fixed_indices] = interpolateNanData(data,varargin)
%
%   
%   seg_worm.feature_helpers.interpolateNanData
%
%
%   Consider using in:
%   seg_worm.feature_helpers.locomotion.getForaging
%   seg_worm.feature_helpers.locomotion.getOmegaAndUpsilonTurns

error('Not yet implemented')

in.good_indices = [];
in.fix_indices  = [];
in.use_extrapolation = false;
in.max_sample_gap = [];
in = sl.in.processVarargin(in,varargin);



end