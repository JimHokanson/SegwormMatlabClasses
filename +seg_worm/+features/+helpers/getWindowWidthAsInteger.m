function window_width_integer = getWindowWidthAsInteger(scale_time,fps)
%
%
%       seg_worm.feature_helpers.getWindowWidthAsInteger
%
%
%   Used by:
%   seg_worm.feature_helpers.path.wormPathCurvature
%   seg_worm.feature_helpers.computeVelocity

% The scale must be odd.
%--------------------------------------------------------------------------
%To this end we take the value obtained and round it down and up, we choose
%the odd value. This is made tricky if it is an even integer to start, such
%that rounding up or down both result in an even integer. In this case we
%add 1
window_width_as_samples = scale_time * fps;

half_scale = round(window_width_as_samples/2);
window_width_integer = 2*half_scale + 1;

%OLD_CODE
%{

scale_low  = floor(window_width_as_samples);
scale_high = ceil(window_width_as_samples);

if scale_low == scale_high
    if mod(scale_low,2) == 0
        window_width_integer = window_width_as_samples + 1;
    else
        window_width_integer = scale_low;
    end
elseif mod(scale_high,2) == 0
    window_width_integer = scale_low;
else
    window_width_integer = scale_high;
end

%}


end