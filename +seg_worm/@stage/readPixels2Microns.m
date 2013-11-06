function readPixels2Microns(obj)
%READPIXELS2MICRONS Read the experiment information file and compute the
%   scale for converting pixels to microns as well as the rotation matrix.
%
%   [PIXEL2MICRONSCALE ROTATION] = READPIXELS2MICRONS(INFOFILE)
%
%   Input:
%       infoFile - the XML file with the experiment information
%
%   Outputs:
%       pixel2MicronScale - the scale for converting pixels to microns
%       rotation          - the rotation matrix
%
% See also PIXELS2MICRONS
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

info = obj.info;
% Compute the microns/pixels.
%pixel2MicronScale = [pixelsX / micronsX, pixelsY / micronsY];
pixel2MicronX     = info.pixels_X / info.microns_X;
pixel2MicronY     = info.pixels_Y / info.microns_Y;
normScale         = sqrt((pixel2MicronX ^ 2 + pixel2MicronY ^ 2) / 2);
obj.pixel_2_micron_scale =  normScale * [sign(pixel2MicronX) sign(pixel2MicronY)];

% Compute the rotation matrix.
%rotation = 1;
angle = atan(pixel2MicronY / pixel2MicronX);
if angle > 0
    angle = pi / 4 - angle;
else
    angle = pi / 4 + angle;
end
cosAngle     = cos(angle);
sinAngle     = sin(angle);
obj.rotation = [cosAngle, -sinAngle; sinAngle, cosAngle];
end
