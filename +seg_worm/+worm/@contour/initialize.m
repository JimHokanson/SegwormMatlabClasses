function initialize(obj,img,frame_number,verbose,num_erode,num_dilate,is_normalized)

%The organization of this code could use a bit of work ...
pe = seg_worm.parse_error(img,frame_number,verbose);

helper__initRoughContour(obj,pe,img,num_erode,num_dilate,is_normalized)
if obj.parse_error, return; end

obj.cleanWorm();

[obj.error_num,obj.error_msg] = pe.contourTooSmall(obj.cleaned_pixels,obj.N_SEGS,img);
obj.parse_error = ~isempty(obj.error_msg);
if obj.parse_error, return; end

helper__computeAngleInfo(obj,pe)

obj.init__getHeadTailI();



keyboard

end

function helper__initRoughContour(obj,pe,img,num_erode,num_dilate,is_normalized)
%
%
%
%   1) Grayscale
%   2) OTSU
%   3) Dilate and erode
%   4) Identification of biggest blob
%   5) Outer contour tracing
%
%   Populates:
%   -------------------------
%   .rough_pixels


% Convert the image to grayscale.
if (size(img,3) == 3)
    img = rgb2gray(img);
end

% Store the original then binarize the image.
img  = seg_worm.cv.otsuImg(img, is_normalized);
img  = ~img; %invert image, normally the pixels are white, the background black
%which is different from

% Erode and dilate the binary image.
%-------------------------------------------
if ~isempty(num_erode) && num_erode > 0
    img = imerode(img, strel('disk', num_erode));
end
if ~isempty(num_dilate) && num_dilate > 0
    img = imdilate(img, strel('disk', num_dilate));
end

obj.image_used = img;

%Find the worm blob
%----------------------------------------
%The largest connected blob is the worm
cc = bwconncomp(img); %img tlbx call
worm_pixels = [];
if ~isempty(cc.PixelIdxList)
    %PixelIdxList - cell array of arrays, each array
    %contains indices of connected pixels
    [~,I] = max(cellfun('length',cc.PixelIdxList));
    worm_pixels = cc.PixelIdxList{I};
end

[obj.error_num,obj.error_msg] = pe.noWormFound(worm_pixels);
obj.parse_error = ~isempty(obj.error_msg);
if obj.parse_error, return; end

% Find a point on the contour.
%----------------------------------------
[y, x] = ind2sub(size(img), min(worm_pixels));
seed   = [x y];

% Trace the contour clockwise.
%----------------------------------------
contour_pixels = seg_worm.cv.bwClockTrace(img,seed,true);

[obj.error_num,obj.error_msg] = pe.contourTouchesBoundary(contour_pixels,img);
obj.parse_error = ~isempty(obj.error_msg);
if obj.parse_error, return; end

[obj.error_num,obj.error_msg] = pe.contourTooSmall(contour_pixels,obj.N_SEGS,img);
obj.parse_error = ~isempty(obj.error_msg);
if obj.parse_error, return; end

obj.rough_pixels = contour_pixels;

end

function helper__computeAngleInfo(obj,pe)

%TODO: Might be nice to smooth in this function ...

import seg_worm.cv.*

cleaned_pixels_local = obj.cleaned_pixels;
% Compute the contour's local high/low-frequency curvature.
% Note: worm body muscles are arranged and innervated as staggered pairs.
% Therefore, 2 segments have one theoretical degree of freedom (i.e. one
% approximation of a hinge). In the head, muscles are innervated
% individually. Therefore, we sample the worm head's curvature at twice the
% frequency of its body.
% Note 2: we ignore Nyquist sampling theorem (sampling at twice the
% frequency) since the worm's cuticle constrains its mobility and practical
% degrees of freedom.
cc_lengths = seg_worm.cv.circComputeChainCodeLengths(cleaned_pixels_local);


%??? - This doesn't make any sense ...
%total_length = cc_lengths(end), so why add the first element?
%
%   %Correct code ...
%   avg_worm_segment_length = cc_lengths(end)/obj.N_SEGS;
%
worm_segment_length = (cc_lengths(1) + cc_lengths(end)) /obj.N_SEGS;

hf_angle_edge_length = worm_segment_length;
hf_angles            = seg_worm.cv.circCurvature(cleaned_pixels_local, hf_angle_edge_length, cc_lengths);

lf_angle_edge_length = 2 * hf_angle_edge_length;
lf_angles            = seg_worm.cv.circCurvature(cleaned_pixels_local, lf_angle_edge_length, cc_lengths);

obj.worm_segment_length = worm_segment_length;
obj.hf_edge_length = hf_angle_edge_length;
obj.lf_edge_length = lf_angle_edge_length;
obj.lf_angles      = lf_angles;
obj.hf_angles_raw  = hf_angles;
obj.cc_lengths     = cc_lengths;

end

