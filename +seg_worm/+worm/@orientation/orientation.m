classdef orientation < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm.orientation
    
    %NOT YET IMPLEMENTED
    %See comments below for content that might go in this class ...
    
    properties
       parent
    end
    
    methods
    end
    
end

%{

% How much confidence do we have in our vulva orientation?
% Note: generally, the vulval side contains less white pixels (a lower
% 50% and 75% CDF for color) and more gray pixels (a lower variance and
% 25% to 75% interquartile range) than the opposing side. We give each
% probability equal weight, then compare. Also, in the absence of
% information, we assume the vulva is on the left side (and use a trick
% to avoid reciprocals in our equations).

isVulvaClockwiseFromHead = 0; % default orientation
vConfidenceScale = 1048576; % 2^20
vConfidence = (rCDF(3) * rCDF(4) * rStdev * (rCDF(4) - rCDF(2))) ...
    / vConfidenceScale;
nvConfidence = (lCDF(3) * lCDF(4) * lStdev * (lCDF(4) - lCDF(2))) ...
    / vConfidenceScale;



headConfidence = struct('head', hConfidence, 'tail', tConfidence);
headOrientation = struct('isFlipped', isHeadTailFlipped, ...
    'confidence', headConfidence);
vulvaConfidence = struct('vulva', vConfidence, 'nonVulva', nvConfidence);
vulvaOrientation = struct('isClockwiseFromHead', isVulvaClockwiseFromHead, ...
    'confidence', vulvaConfidence);
orientation = struct('head', headOrientation, 'vulva', vulvaOrientation);





%}
