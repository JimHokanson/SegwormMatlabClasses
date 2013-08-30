function init__getHeadTailI(obj)

helper__processPossibleAngleErrors(obj,lf_max_peaks,hf_max_peaks)

helper__getHeadTailMain(obj,lf_max_peaks,hf_max_peaks,lf_max_peaks_I,hf_max_peaks_I)

% Orient the contour and angles at the maximum curvature (the head or tail).
if headI > 1
    %IMPORTANT: This redefines properties
    helper__headTailPropSwitch(obj);
end



end

function [lf_max_peaks,hf_max_peaks,lf_max_peaks_I,hf_max_peaks_I] = helper__getLocalPeaks(obj)
%
%
%We'll filter the high frequency data using a moving-average filter to get
%better high frequency estimates. The high

cleaned_pixels_local = obj.cleaned_pixels;

%NOTE: Consider making this a property of the class
worm_seg_size      = size(cleaned_pixels_local,1)/obj.N_SEGS;
hf_angle_edge_size = worm_seg_size;


hf_blur_size       = ceil(hf_angle_edge_size/2);
hf_blur_win        = 1/hf_blur_size*ones(1,hf_blur_size);
hf_angles_smoothed = circConv(obj.hf_angles, hf_blur_win);

%Compute the contour's local high/low-frequency curvature maxima.
%Improvement, could set thresholds for peaks in function
%i.e. only return peaks within a given range ...
[hf_max_peaks,hf_max_peaks_I] = maxPeaksCircDist(hf_angles_smoothed, obj.hf_edge_length, obj.cc_lengths);
[lf_max_peaks,lf_max_peaks_I] = maxPeaksCircDist(obj.lf_angles, obj.lf_edge_length, obj.cc_lengths);


end

function helper__getHeadTailMain(obj)
%Identification of head and tail
%--------------------------------------------------------------------------
if lf_HT_size > 1
    
    %Should be two based on previous error checking
    
    
    % Find the head and tail convexities in the low-frequency sampling.
    % Note: the tail should have a sharper angle.
    
    lf_HT_peaks = lf_max_peaks(lf_HT_I);
    if lf_HT_peaks(1) <= lf_HT_peaks(2)
        head_I = lf_HT_I(1);
        tail_I = lf_HT_I(2);
    else
        head_I = lf_HT_I(2);
        tail_I = lf_HT_I(1);
    end
    
    % Localize the head by finding its nearest, sharpest (but blurred),
    % high-frequency convexity.
    
    dhfHeadI    = abs(cc_lengths_local(head_I) - cc_lengths_local(hf_HT_I));
    dhfHeadI    = min(dhfHeadI, cc_lengths_local(end) - dhfHeadI);
    [~, temp_I] = min(dhfHeadI);
    head_I      = hf_HT_I(temp_I);
    
    % Localize the tail by finding its nearest, sharpest (but blurred),
    % high-frequency convexity.
    dhfTailI    = abs(cc_lengths_local(tail_I) - cc_lengths_local(hf_HT_I));
    dhfTailI    = min(dhfTailI, cc_lengths_local(end) - dhfTailI);
    [~, temp_I] = min(dhfTailI);
    tail_I      = hf_HT_I(temp_I);
    
    % The high-frequency sampling identifies the head and tail.
elseif hf_HT_size < 3
    
    % Find the head and tail convexities in the high-frequency sampling.
    % Note: the tail should have a sharper angle.
    mhfHTI = hf_max_peaks_I(hf_HT_I);
    mhfHTP = hf_max_peaks(hf_HT_I);
    if mhfHTP(1) <= mhfHTP(2)
        head_I = mhfHTI(1);
        tail_I = mhfHTI(2);
    else
        head_I = mhfHTI(2);
        tail_I = mhfHTI(1);
    end
    
    % The high-frequency sampling identifies several, potential heads/tails.
else
    
    % Initialize our head and tail choicse.
    mhfHTI  = hf_max_peaks_I(hf_HT_I);
    mhfHTI1 = mhfHTI(1);
    mhfHTI2 = mhfHTI(2);
    
    % How far apart are the head and tail?
    dmhfHTI12 = abs(cc_lengths_local(mhfHTI(1)) - cc_lengths_local(mhfHTI(2)));
    dmhfHTI12 = min(dmhfHTI12, cc_lengths_local(end) - dmhfHTI12);
    
    % Search for the 2 sharp convexities that are furthest apart.
    for i = 1:(hf_HT_size - 1)
        for j = (i + 1):hf_HT_size
            
            % How far apart are these 2 convexities?
            dmhfHTIij = abs(cc_lengths_local(mhfHTI(i)) - ...
                cc_lengths_local(mhfHTI(j)));
            dmhfHTIij = min(dmhfHTIij, cc_lengths_local(end) - dmhfHTIij);
            
            % These 2 convexities are better head and tail choices.
            if dmhfHTIij > dmhfHTI12
                mhfHTI1 = mhfHTI(i);
                mhfHTI2 = mhfHTI(j);
                dmhfHTI12 = dmhfHTIij;
            end
        end
    end
    
    % Which convexity is the head and which is the tail?
    % Note: the tail should have a sharper angle.
    if mhfCAngles(mhfHTI1) < mhfCAngles(mhfHTI2)
        head_I = mhfHTI1;
        tail_I = mhfHTI2;
    else
        head_I = mhfHTI2;
        tail_I = mhfHTI1;
    end
end

cc_lengths_local = obj.cc_lengths;

obj.head_I = head_I;
obj.tail_I = tail_I;

% Find the length of each side.
if head > tailI
    size_1 = cc_lengths_local(headI) - cc_lengths_local(tailI);
    size_2 = cc_lengths_local(end)   - cc_lengths_local(headI) + cc_lengths_local(tailI);
else
    size_1 = cc_lengths_local(tailI) - cc_lengths_local(headI);
    size_2 = cc_lengths_local(end)   - cc_lengths_local(tailI) + cc_lengths_local(headI);
end

% Are the sides within 50% of each others size?
% Note: if a worm's length from head to tail is at least twice larger
% on one side (relative to the other), than the worm must be touching
% itself.
if min(size_1, size_2)/ max(size_1, size_2) <= .5
    obj.parse_error = true;
    obj.error_num   = 106;
    obj.error_msg   = ['The worm length, from head to tail, is more than ' ...
        'twice as large on one side than it is on the other. ' ...
        'Therefore, the worm is coiled or obscured and cannot be segmented.'];
    return;
end

end

function helper__headTailPropSwitch(obj)
%
%
%   This occurs if head_I is not the first pixel
%


%????? - what is this reallly doing????

cleaned_pixels_local = obj.cleaned_pixels;
cc_lengths_local     = obj.cc_lengths;
hf_angles_local      = obj.hf_angles_raw;
lf_angles_local      = obj.lf_angles;
head_I_local         = obj.head_I;
tail_I_local         = obj.tail_I;

cleaned_pixels_local = [...
    cleaned_pixels_local(head_I_local:end,:); ...
    cleaned_pixels_local(1:(head_I_local - 1),:)];
cc_lengths_local = [...
    cc_lengths_local(headI:end) - cc_lengths_local(head_I_local - 1); ...
    cc_lengths_local(1:(headI - 1)) + ...
    (cc_lengths_local(end) - cc_lengths_local(head_I_local - 1))];
hf_angles_local = [hf_angles_local(head_I_local:end); hf_angles_local(1:(head_I_local - 1))];
lf_angles_local = [lf_angles_local(head_I_local:end); lf_angles_local(1:(head_I_local - 1))];

%What is this doing????
lfCMaxI = lfCMaxI - head_I_local + 1;
wrap = lfCMaxI < 1;
lfCMaxI(wrap) = lfCMaxI(wrap) + length(lfCAngles);
tail_I_local = tail_I_local - head_I_local + 1;
head_I_local = 1;
if tail_I_local < 1
    tail_I_local = tail_I_local + size(contour_obj, 1);
end

obj.clean_pixels = cleaned_pixels_local;
obj.cc_lengths   = cc_lengths_local;
obj.hf_angles    = hf_angles_local;
obj.lf_angles    = lf_angles_local;
obj.head_I       = head_I_local;
obj.tail_I       = tail_I_local;

end

function helper__processPossibleAngleErrors(obj,lf_max_peaks,hf_max_peaks)

LF_ANGLE_CUTOFF = 90;
HF_ANGLE_CUTOFF = 60;

%ERRORS
%--------------------------------------------------------------------------
% Are there too many possible head/tail points?
lf_HT_I    = lf_max_peaks_I(lf_max_peaks > LF_ANGLE_CUTOFF); %lf_HT - low frequency, head/tail
lf_HT_size = length(lf_HT_I);
if lf_HT_size > 2
    obj.parse_error = true;
    obj.error_num = 104;
    obj.error_msg = ['The worm has 3 or more low-frequency sampled convexities' ...
        'sharper than 90 degrees (possible head/tail points).'];
    return
end

hf_HT_I    = hf_max_peaks_I(hf_max_peaks > HF_ANGLE_CUTOFF);
hf_HT_size = length(hf_HT_I);
if hf_HT_size < 2
    obj.parse_error = true;
    obj.error_num = 105;
    obj.error_msg = ['The worm contour has less than 2 high-frequency sampled '...
        'convexities sharper than 60 degrees (the head and tail). ' ...
        'Therefore, the worm is coiled or obscured and cannot be segmented.'];
    return
end


end