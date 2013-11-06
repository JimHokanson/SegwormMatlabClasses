function initTailHeadIndices(obj)
%
%
%   seg_worm.worm.contour.initTailHeadIndices

mask_hf   = obj.hf_ap_max > obj.HF_MAX_ANGLE_LIMIT;
mask_lf   = obj.lf_ap_max > obj.LF_MAX_ANGLE_LIMIT;

hf_HT_I   = obj.hf_ap_max_I(mask_hf);
lf_HT_I   = obj.lf_ap_max_I(mask_lf);

n_hf_HT_I = length(hf_HT_I);
n_lf_HT_I = length(lf_HT_I);

hf_HT_max_peaks = obj.hf_ap_max(mask_hf);
lf_HT_max_peaks = obj.lf_ap_max(mask_lf);

%ERRORS
%--------------------------------------------------------------------------
% Are there too many possible head/tail points?
obj.error_handler.insufficientHTOptions(n_lf_HT_I,n_hf_HT_I);
if obj.parse_error, return; end

%--------------------------------------------------------------------------
if n_lf_HT_I > 1
    [head_I,tail_I] = helper__getHeadTailByLF(obj,lf_HT_max_peaks,lf_HT_I,hf_HT_I);
elseif n_hf_HT_I < 3
    % The high-frequency sampling identifies the head and tail.
    [head_I,tail_I] = helper__getHeadTailByHF(hf_HT_max_peaks,hf_HT_I);
else
    % The high-frequency sampling identifies several, potential heads/tails.
    [head_I,tail_I] = helper__getHeadTailByManyHF;
end

obj.head_I = head_I;
obj.tail_I = tail_I;

obj.error_handler.lopsidedSides(obj);
if obj.parse_error, return; end

% Orient the contour and angles at the maximum curvature (the head or tail).
if head_I > 1
    %???? Why would head_I be at 1, what about our contour parsing makes
    %this likely
    %IMPORTANT: This redefines properties
    helper__headTailPropSwitch(obj);
end



end

function [head_I,tail_I] = helper__getHeadTailByLF(obj,lf_HT_max_peaks,lf_HT_I,hf_HT_I)

    cc_lengths_local = obj.cc_lengths;

    % Find the head and tail convexities in the low-frequency sampling.
    % Note: the tail should have a sharper angle.
    
    if lf_HT_max_peaks(1) <= lf_HT_max_peaks(2)
        lf_head_I = lf_HT_I(1);
        lf_tail_I = lf_HT_I(2);
    else
        lf_head_I = lf_HT_I(2);
        lf_tail_I = lf_HT_I(1);
    end
    
    % Localize the head by finding its nearest, sharpest (but blurred),
    % high-frequency convexity.
    %
    %i.e.
    %Find the closest hf_HT_I to lf_head_I
    %NOTE: The distance could be closest by "left (wrap around)" or "right"

    %NOTE: The min is over the array of options     
    [~, temp_I] = min(helper_getMinCCDistance(cc_lengths_local,lf_head_I,hf_HT_I));
    head_I      = hf_HT_I(temp_I);
    
    % Localize the tail by finding its nearest, sharpest (but blurred),
    % high-frequency convexity.
    [~, temp_I] = min(helper_getMinCCDistance(cc_lengths_local,lf_tail_I,hf_HT_I));
    tail_I      = hf_HT_I(temp_I);
end

function [head_I,tail_I] = helper__getHeadTailByHF(hf_HT_max_peaks,hf_HT_I)

    % Find the head and tail convexities in the high-frequency sampling.
    % Note: the tail should have a sharper angle.
    if hf_HT_max_peaks(1) <= hf_HT_max_peaks(2)
        head_I = hf_HT_I(1);
        tail_I = hf_HT_I(2);
    else
        head_I = hf_HT_I(2);
        tail_I = hf_HT_I(1);
    end

end

function [head_I,tail_I] = helper__getHeadTailByManyHF(obj,hf_HT_max_peaks,hf_HT_I)

    cc_lengths_local = obj.cc_lengths;

    n_peaks = length(hf_HT_max_peaks);
    combos  = nchoosek(1:n_peaks,2);
    
    dist_all = arrayfun(@(x,y) helper_getMinCCDistance(cc_lengths_local,...
        hf_HT_I(combos(:,1)),hf_HT_I(combos(:,2))));
    
    %This isn't finished, but is close
    keyboard
    
% % %     % Initialize our head and tail choicse.
% % %     hf_HT_I1 = hf_HT_I(1);
% % %     hf_HT_I2 = hf_HT_I(2);
% % %     
% % %     % How far apart are the head and tail?
% % %     dmhfHTI12 = abs(cc_lengths_local(hf_HT_I(1)) - cc_lengths_local(hf_HT_I(2)));
% % %     dmhfHTI12 = min(dmhfHTI12, cc_lengths_local(end) - dmhfHTI12);
% % %     
% % %     % Search for the 2 sharp convexities that are furthest apart.
% % %     for i = 1:(hf_HT_size - 1)
% % %         for j = (i + 1):hf_HT_size
% % %             
% % %             % How far apart are these 2 convexities?
% % %             dmhfHTIij = abs(cc_lengths_local(hf_HT_I(i)) - ...
% % %                 cc_lengths_local(hf_HT_I(j)));
% % %             dmhfHTIij = min(dmhfHTIij, cc_lengths_local(end) - dmhfHTIij);
% % %             
% % %             % These 2 convexities are better head and tail choices.
% % %             if dmhfHTIij > dmhfHTI12
% % %                 hf_HT_I1 = hf_HT_I(i);
% % %                 hf_HT_I2 = hf_HT_I(j);
% % %                 dmhfHTI12 = dmhfHTIij;
% % %             end
% % %         end
% % %     end
    
    % Which convexity is the head and which is the tail?
    % Note: the tail should have a sharper angle.
    if hf_HT_max_peaks(hf_HT_I1) < hf_HT_max_peaks(hf_HT_I2)
        head_I = hf_HT_I1;
        tail_I = hf_HT_I2;
    else
        head_I = hf_HT_I2;
        tail_I = hf_HT_I1;
    end
end

function helper__headTailPropSwitch(obj)
%
%
%   This occurs if head_I is not the first pixel

head_I_local = obj.head_I;

cur_prop_value = obj.pixels;
obj.pixels = [cur_prop_value(head_I_local:end,:); cur_prop_value(1:(head_I_local - 1),:)];

obj.computeAngleInfo();

tail_I_local = obj.tail_I;
tail_I_local = tail_I_local - head_I_local + 1;
if tail_I_local < 1
   tail_I_local = tail_I_local + obj.n_pixels;
end
obj.head_I = 1;
obj.tail_I = tail_I_local;

end

function min_dist = helper_getMinCCDistance(cc_lengths_local,index1,index_or_indices)
%
%   min_dist = helper_getMinCCDistance(cc_lengths_local,index1,index_or_indices)
%
%   NOTE: index1 and index_or_indices both index into cc_lengths_local
%
%   OUTPUT:
%   min_dist: [1 x n], for each index rea
%
%   NOTE: I should should check one input is singular, both can't be
%   arrays. In reality it doesn't matter which is the array

    %Get's the distance between two points, taking into account the
    %circular nature of the location and that you could go one way or the
    %other to get from one point to the other
    dRight   = abs(cc_lengths_local(index1) - cc_lengths_local(index_or_indices));
    dLeft    = cc_lengths_local(end) - dRight;
    min_dist = min(dRight,dLeft);
end
