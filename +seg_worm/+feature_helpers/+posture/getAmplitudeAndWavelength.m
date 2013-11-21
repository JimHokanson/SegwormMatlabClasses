function [amplitude,wavelength,trackLength] = getAmplitudeAndWavelength(theta,xx,yy,wormLen)
%
%   JAH NOTE: Schafer lab can't find this code, I might need to reimplement
%   it ...
%   
%   seg_worm.feature_helpers.posture.getAmplitudeAndWavelength
%
%   amplitude
%
% worm.posture.amplitude.max        = ampt;
% worm.posture.amplitude.ratio      = ampRatio;
% worm.posture.wavelength.primary   = wavelnths(1);
% worm.posture.wavelength.secondary = wavelnths(2);
% worm.posture.tracklength = trackLen;

%Note from old code
%--------------------------------------------------------------------------
% Note: This code was taken directly from matlab code published by Cronin
% et al. in this publication:
% BMC Genetics, 2005
% C.J. Cronin, J.E. Mendel, S. Mukhtar, Young-Mee Kim, R.C. Stirb, J. Bruck,
% P.W. Sternberg
% "An automated system for measuring parameters of nematode
% sinusoidal movement" BMC Genetics 2005, 6:5
%
% It was received from Cronin from personal communication and his agreement
% to use it was received.
% It is therefore not owend or claimed to be owned by creators of the Worm
% Analysis Toolbox at Schafer lab, MRC Laboratory of Molecular Biology


%JAH: This code is a work in progress ...

%TODO: Redo code for all frames ...


%----TRACK AMPLITUDE--------------------------------------

% Define the rotation matrix, 
% The direction of vector rotation is counterclockwise if theta is positive 
% (e.g. 90°), and clockwise if theta is negative (e.g. -90°).
B = [cos(-theta)  -sin(-theta)   ;
     sin(-theta)   cos(-theta)   ];
 
% Rotate the skeleton coordinates
ww = B*[xx';yy'];

% Define number of points and intervals in the skeleton
npts       = size(ww, 2);
nintervals = npts - 1;

% Parse transformed worm-coordinate matrix 
% x&y coordinate vectors
wwx = ww(1,:);
wwy = ww(2,:);

%Amplitude Calculations
%--------------------------------------------------------------------------
% We have rotated the coordinates around the origin at the axes, not at the
% centroid of the worm. Now center at the origin by subtracting the
% centroid coordinates

wwx = wwx - nanmean(wwx); 
wwy = wwy - nanmean(wwy);
% Calculate track amplitude 
amplitude.max = nanmax(wwy) - nanmin(wwy);

% Calculate track length
trackLength = nanmax(wwx)-nanmin(wwx);
amp1 = max(wwy(wwy>0));
amp2 = max(abs(wwy(wwy<0)));
amplitude.ratio = min(amp1,amp2)/max(amp1,amp2);
%worm.posture.amplitude.ratio





%Correction Orientaton?
%------------------------------------------------------------------






















%MISSING FILE
%--------------------------------------------------------------------------
%getAmpWavelength
%
% BMC Genetics, 2005
% C.J. Cronin, J.E. Mendel, S. Mukhtar, Young-Mee Kim, R.C. Stirb, J. Bruck,
% P.W. Sternberg
% "An automated system for measuring parameters of nematode
% sinusoidal movement" BMC Genetics 2005, 6:5
%
%   Code supposedly at:
%   http://wormlab.caltech.edu/publications/CaltechTracker(BMCGenetics2005).zip
%
%   Doesn't contain getAmpWavelength
%--------------------------------------------------------------------------



%From Sternberg code, might use:
%tracks.m           - ampt, wavelnth
%
%   ??? What is trackLen and ampRatio

%   -> ampt - maximum amplitude along the worm body
%   -> ampRatio - maximum on opposiding sides, smaller is numerator
%   -> uses orientation from ellipse
%   -> trackLen - horizontal length of projection
%       -> I think this is the opposite of ampt
%
%       trackLen/ampt -> ampRatio???

%wavelength - x must be monotonic after rotation
%primary - peak of fft
%secondary - 2nd peak of fft, must exceed 0.5 * 1st peak
%
%   NOTE: Wavelength can never exceed 2x worm length
%   -> If greater then what? NaN?
%
%
%   ???? - velocity based, but their code isn't ...

%MISSING FUNCTION
[ampt, wavelnths, trackLen, ampRatio] = getAmpWavelength(thetaValRad, skCoords, wormLen, mydata.mainImg, guiItem, 0, 0);

worm.posture.amplitude.max        = ampt;
worm.posture.amplitude.ratio      = ampRatio;
worm.posture.wavelength.primary   = wavelnths(1);
worm.posture.wavelength.secondary = wavelnths(2);
worm.posture.tracklength = trackLen;