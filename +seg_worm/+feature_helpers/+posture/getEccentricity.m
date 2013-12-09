function [eccentricity, orientation] = getEccentricity(xOutline, yOutline, gridSize)
% getEccentricity   
%
%   [eccentricity, orientation] = seg_worm.feature_helpers.posture.getEccentricity(xOutline, yOutline, gridSize)
%
%   Given x and y coordinates of the outline of a region of interest, fill
%   the outline with a grid of evenly spaced points and use these to
%   calculate the eccentricity and orientation of the equivalent ellipse.
%
%   Placing points in the contour is a well known computer science problem
%   known as the Point-in-Polygon problem.
%
%   http://en.wikipedia.org/wiki/Point_in_polygon
%
%   This function became a lot more complicated in an attempt to make it 
%   go much faster. The complication comes from the simplication that can
%   be made when the worm doesn't bend back on itself at all.
%
%
%   OldName: getEccentricity.m
%
%   Inputs:
%   =======================================================================
%   xOutline : The x coordinates of the contour. In particular the contour
%               starts at the head and goes to the tail and then back to
%               the head (although no points are redundant)
%   yOutline : The y coordinates of the contour "  "
%   gridSize : The # of points to place in the long dimension. More points
%              gives a more accurate estimate of the ellipse but increases
%              the calculation time.
%
%   Outputs:
%   =======================================================================
%   eccentricity - The eccentricity of the equivalent ellipse
%   orientation  - The orientation angle of the equivalent ellipse
%
%   Nature Methods Description
%   =======================================================================
%   Eccentricity. 
%   ------------------
%   The eccentricity of the worm’s posture is measured using
%   the eccentricity of an equivalent ellipse to the worm’s filled contour.
%   The orientation of the major axis for the equivalent ellipse is used in
%   computing the amplitude, wavelength, and track length (described
%   below).


MIDDLE_INDEX = 49;

%center the worm at its centroid
%--------------------------------------------------------------------------
xOutline_mc = bsxfun(@minus,xOutline,mean(xOutline,1)); %_mc - mean centered
yOutline_mc = bsxfun(@minus,yOutline,mean(yOutline,1));


%Rough rotation of the worm
%--------------------------------------------------------------------------
%
%   We rotate the worm by the vector formed from going from the first point
%   to the last point. This accomplishes two things.
%
%   1) For the brute force approach, this makes it so the worm is
%   encompassed better by a rectangle instead of a square, which means that
%   there are fewer grid points to test that are not in the worm (as
%   opposed to a worm that is on a 45 degree angle that would have many
%   points on the grid outside of the worm). In other words, the bounding
%   box is a better approximation of the worm when the worm is rotated then
%   when it is not. In the brute force case we place points in the bounding
%   box, so the smaller the bounding box, the faster the code will run.
%
%   2) This allows us to hardcode only looking for "simple" worms 
%   (see description below) that vary in the x-direction. Otherwise we
%   might need to also look for simple worms in the y direction. Along
%   similar lines this makes it more likely that we will get simple worms.

y = yOutline_mc(MIDDLE_INDEX,:) - yOutline_mc(1,:); %_mc - mean centered
x = xOutline_mc(MIDDLE_INDEX,:) - xOutline_mc(1,:);

rot_angle = atan2(y,x);

%I expanded the rotation matrix to allow processing all frames at once
%
%   i.e. rather than R*X (Matrix multiplication)
%
%   I have r(1)*X(1) + r(2)*X(2) etc, but where X(1), X(2), etc is really
%   a vector of values, not just a singular value like you would need for
%   R*X
%
%_rot - rotated
xOutline_mc_rot = bsxfun(@times,xOutline_mc,cos(rot_angle))  + bsxfun(@times,yOutline_mc,sin(rot_angle));
yOutline_mc_rot = bsxfun(@times,xOutline_mc,-sin(rot_angle)) + bsxfun(@times,yOutline_mc,cos(rot_angle));

%TODO: replace refs of _mc with _rot, this is a temporary fix
xOutline_mc = xOutline_mc_rot; 
yOutline_mc = yOutline_mc_rot;

%{
%Debugging code ...
%Examine the rotated worms 
%-------------------------------------------------
for iFrame = 1:10:size(xOutline_mc,2)
   scatter(xOutline_rot(:,iFrame),yOutline_rot(:,iFrame))
   title(sprintf('Frame %d',iFrame));
   axis equal
   pause
end

%}

%Determine 'simple' worms
%--------------------------------------------------------------------------
%
%   Simple worms have two sides that go from - to + (or + to -), with no bending
%   backwards. 
%
%   e.g.
%
%   x_side_1 = 1,2,3 - simple
%   x_side_1 = 1 3 2 - not simple
%
%   If we find a simple worm then we can interpolate the y
%   values at the x-grid locations. This removes the need to sort the x
%   values in the grid with respect to the x-values of the contour. Once
%   we have the x-values aligned, we round the y-values to grid points.
%   Doing this allows us to simply count the points off that are between
%   the two y-values.
%
%   i.e., let's say we have x_side_1 = 0, 1.2, 2.3, 3.5
%
%   Normally a grid algorithm would first have to find which set of x side
%   values the grid point is in between. A naive point in polygon algorithm
%   will not know that a series of grid points with the same x values and
%   differing y-values all have the same x bounds (from the polygon or in
%   this case contour).
%
%   so, instead we can say, well our x grid values are at:
%   0,1,2,3,4 - note the event spacing
%
%   and then for each of those grid values, we can get the y-values of the
%   contour
%
%   For the y-values (interpolated contour points), consider the points 13
%   and 33, and grid values at every 5, starting at 11 (minimum extent of
%   the contour). We take 13 and round up to the nearest grid point, which
%   would be 16. We take 33 and round down to the nearest grid point, which
%   would be floor((33 - 11)/5) => 20, 20 + 11 => 31. So now that we have
%   16 and 31 the points between 13 and 33 are: 16:5:31 or 16, 21, 26, 31
%
%   This is illustrated below:
%
%           33 - y value of higher contour
%       
%           13 - y value of lower contour
%           11 - min
%
%
%   If 11 is min, and 5 is our spacing, then our grid points to test will
%   be 11,16,21,26,31,36,41, etc
%
%   If we round 13 up and 33 down to their appropriate locations (16 and
%   31), then we know that at this x value, the grid points 16:5:31 will
%   all be in the worm

TOL = 1.0e-12; %We need a bit of a tolerance on whether or not a point is 
%in the worm or outside the worm

%X extents of worms
[mx_x,mx_x_I] = max(xOutline_mc_rot,[],1);
[mn_x,mn_x_I] = min(xOutline_mc_rot,[],1);

%We don't know if the worm will be going from left to right or right to
%left, we need slightly different code later on depending on which
%
%NOTE: we are testing the indices, not the values
min_first_mask = mn_x_I < mx_x_I;
min_last_mask  = mx_x_I < mn_x_I;

%Y extents
mx_y = max(yOutline_mc,[],1);
mn_y = min(yOutline_mc,[],1);

%We take the difference between worm contour points in the x dimension
%Since the contour is circular, we compute the difference between the 
%first and last points
%
%NOTE: 2 - 1, 3 - 2, ..., N - (N-1), then we need 1 - N
d = [diff(xOutline_mc,1,1); xOutline_mc(1,:) - xOutline_mc(end,:)];

n_frames = size(xOutline,2);
n_c      = size(xOutline,1);

y_interp_1 = NaN(gridSize,n_frames); %Tentatively the bottom of the worm, we'll
%explicitly check this out later ...
y_interp_2 = NaN(gridSize,n_frames);
x_interp   = NaN(gridSize,n_frames);

%--------------------------------------------------------------------------
%NOTE: The basic difference between these loops is in how x1 and x2 are
%defined. For interpolation we must always go from a lower value to a
%higher value (to use the quick method of interpolation in Matlab). Note, 
%the sign on the comparison is also different.

%In these two loops we calculate the interpolated y contour values at the
%x-grid locations
%x1: from negative to postive x
%x2: from positive to negative x - needs to be reversed during interpolation
for iFrame = find(min_first_mask)
    
    %NOTE: x1 and x2 are both inclusive of the min and max to allow
    %the interpolation code to work without extrapolation
    x1 = mn_x_I(iFrame):(mx_x_I(iFrame));
    
    if all(d(x1(1:end-1),iFrame) > 0)
        
        x2 = [mx_x_I(iFrame):n_c    1:(mn_x_I(iFrame))];
        
        if all(d(x2(1:end-1),iFrame) < 0)
            [y_interp_1(:,iFrame),y_interp_2(:,iFrame),x_interp(:,iFrame)] = ...
                h__getInterpValues(x1,x2,xOutline_mc,yOutline_mc,gridSize,iFrame,mn_x_I,mx_x_I);
        end
    end
end

%--------------------------------------------------------------------------
%x1 and x2 maintain their definitions, but they are calculated differently
for iFrame = find(min_last_mask)
    
    %NOTE: min last means that the index of the maximum x value is less
    %than the index of the minimum x value, i.e. we might have values
    %such as mx_x_I = 10  and   mn_x_I = 50
    %
    %but note that the x values themselves might be going from 100 to 1
    x2 = mx_x_I(iFrame):(mn_x_I(iFrame));
    
    if all(d(x2(1:end-1),iFrame) < 0)
        
        x1 = [mn_x_I(iFrame):n_c    1:(mx_x_I(iFrame))];
        
        if all(d(x1(1:end-1),iFrame) > 0)
            [y_interp_1(:,iFrame),y_interp_2(:,iFrame),x_interp(:,iFrame)] = ...
                h__getInterpValues(x1,x2,xOutline_mc,yOutline_mc,gridSize,iFrame,mn_x_I,mx_x_I);
        end
    end
end

%Computing the y-values of the grid points
%--------------------------------------------------------------------------
%   This bit is a little tricky. We want to place our grid locations
%   at each dy, using the appropriate offset from the minimum of the worm.
%
%   1) - subtract minimum for centering
%   2) - round y1 up to nearest dy and y2 down to nearest dy, this is done
%   by  y1 = ceil(y1/dy)*dy
%       y2 = floor(y1/dy)*dy
%       
%       These lines are somewhat messy because in Matlab we tradeoff
%       efficiency of calculations for multiple lines. Since these are
%       relatively slow and the memory allocation can take a decent amount
%       of time I've tried to condense lines rather than having one
%       calculation per line. We also get a bit awkward because of the
%       bsxfun notation.
%
%   3) - read minimum
%
%   The minimum is subtracted to allow better coverage of the worm. We
%   don't want do test each dy 0, dy, 2dy, instead we want it so that the
%   extreme y-value of the worm is the starting point
%
%   min_y, min_y + dy, min_y + 2dy, etc
%
%   This is why we must first subtract the minimum. We then readd it for 
%   doing centroid calculations, although I think this step may be
%   unnecessary ...


%A simple worm is one in which after rotation the x position is always
%increasing or decreasing. If this is the case, then we will have defined
%the new 'x_interp' variable, if not, it will still be NaN. In Matlab, any
%is not true for NaN, for all it is ...
simple_worm_mask = any(x_interp,1);

dy = x_interp(2,:) - x_interp(1,:); %Same as dx to have "grid" sampling
min_y_interp = min(y_interp_1,[],1); 

%Filter values to only do calculations on simple worm frames
dy = dy(simple_worm_mask);
min_y_interp = min_y_interp(simple_worm_mask);

%NOTE: The tolerance was subtracted to try and ensure points exactly
%on the borders of the worm are included as well
%
%NOTE: The rounded values are only a subset of the original frames, so we
%need to be careful about which index we are using
y_1_rounded = bsxfun(@minus,y_interp_1(:,simple_worm_mask),min_y_interp)-TOL;
y_2_rounded = bsxfun(@minus,y_interp_2(:,simple_worm_mask),min_y_interp)-TOL;
y_1_rounded = bsxfun(@times,ceil(bsxfun(@rdivide,y_1_rounded,dy)),dy);
y_2_rounded = bsxfun(@times,floor(bsxfun(@rdivide,y_2_rounded,dy)),dy);
y_1_rounded = bsxfun(@plus,y_1_rounded,min_y_interp);
y_2_rounded = bsxfun(@plus,y_2_rounded,min_y_interp);

%size of y_1_rounded and y_2_rounded - [gridSize x n_frames]

%??? Does aligning to the minimum bias the result??? We should get an
%equally valid answer if we just round to zero ...
%
%JAH: I might change this to remove the subtraction of the minimums

keyboard

%Ensure that y1 < y2 for all frames (DONE), if not then flip (NYI)
%--------------------------------------------------------------------------
%NOTE: we skip the first and last point because they should be the same, although after rounding (floor, ceil)
%they actually will tend to be apart by the amount 'dy', with the opposite relationship that the rest of the data has
%I also filter this out in the loop below by skipping the
%first and last grid point
%
%We need to ensure y1 is less than y2, because below we create vectors
%going from low to high, and if these are reversed, then the vector will
%be empty and the result empty
%
%NOTE: For simplicity I do this calculation on the original data so that I
%don't have to deal with the rounding, which could flip close values (a
%flip indicates that no y grid points are within that x-grid point location
%being tested)
is_y1_lower = all(y_interp_1(2:end-1,simple_worm_mask) < y_interp_2(2:end-1,simple_worm_mask),1);

%NOTE: This can be fixed easily. Any violating frames just need to have y1
%and y2 swapped
if ~all(is_y1_lower)
   error('Assumption violated, y1 is not always less than y2') 
end

%Initialize outputs of the loop
%--------------------------------------------------------------------------
eccentricity = NaN(1,n_frames);
orientation  = NaN(1,n_frames);
%These are temporary arrays for holding the location of grid points that
%fit inside the worm. They are a linerization of all points, so they don't
%have a second dimension, we just pile new points from a worm frame onto
%any old points from that frame.
x_all = zeros(1,gridSize*gridSize);
y_all = zeros(1,gridSize*gridSize);

cur_simple_I = 0; %We have two indices, the original frame numbers and the filtered frame numbers
for iFrame = find(simple_worm_mask) %Original index
    count = 0;
    cur_simple_I = cur_simple_I + 1; %Filtered index
    
    cur_d_unit = dy(cur_simple_I);
    
    %For each x position, we increment from the minimum y value at that x location
    %to the maximum at that location, in the specified steps. We need
    %to hold onto the values for doing the eccentricity and orientation
    %calculations.
    %
    %NOTE: First and last grid points will not contain useful data
    for iIndex = 2:gridSize-1
        
        %Generate appropriate y-values on grid
        temp = y_1_rounded(iIndex,cur_simple_I):cur_d_unit:y_2_rounded(iIndex,cur_simple_I);
        
        %and store ...
        y_all(count+1:count+length(temp)) = temp;
        x_all(count+1:count+length(temp)) = x_interp(iIndex,iFrame);
        count = count + length(temp);
    end

    [eccentricity(iFrame),orientation(iFrame)] = h__calculateSingleValues(x_all(1:count),y_all(1:count));    
end

%{

%Debugging code

plot(xOutline_mc(:,iFrame),yOutline_mc(:,iFrame),'k')
hold on
scatter(x_out_all(:,iFrame),y_interp_1(:,iFrame),'g')
scatter(x_out_all(:,iFrame),y_interp_2(:,iFrame)+0.1,'r')

scatter(x_out_all(:,iFrame),y_1_rounded(:,iFrame),'m')
scatter(x_out_all(:,iFrame),y_2_rounded(:,iFrame)+0.1,'m')

scatter(x_all(1:count),y_all(1:count),'b')
hold off
axis equal

%}

%Use slow grid method for all unfinished worms
%--------------------------------------------------------------------------
xRange_all = mx_x - mn_x;
yRange_all = mx_y - mn_y;
gridAspectRatio_all = xRange_all./yRange_all;

run_mask = ~simple_worm_mask & ~isnan(gridAspectRatio_all);

[eccentricity,orientation] = h__getEccentricityAndOrientation(...
    xOutline_mc,yOutline_mc,xRange_all,yRange_all,gridAspectRatio_all,gridSize,eccentricity,orientation,run_mask);

%Fix the orientation ...
%--------------------------------------------------------------------------
orientation3 = orientation + rot_angle*180/pi;
orientation3(orientation3 > 90)  = orientation3(orientation3 > 90) - 180;
orientation3(orientation3 < -90) = orientation3(orientation3 < -90) + 180;

orientation = orientation3;

end

function [y_interp_1,y_interp_2,x_out_all] = h__getInterpValues(x1,x2,xOutline_mc,yOutline_mc,gridSize,iFrame,mn_x_I,mx_x_I)
%
%
%   [y_interp_1,y_interp_2,x_out_all] = h__getInterpValues(x1,x2,xOutline_mc,yOutline_mc,gridSize,iFrame,mn_x_I,mx_x_I)
%
%
%   This function computes the interpolated y-values on the contour
%   for the x grid locations that we will be testing at. It also computes
%   the x grid (TODO: Might move this ...)
%   
%   Inputs
%   =======================================================================
%   x1          : [1 x m] contour indices with values going from low to high
%   x2          : [1 x n] contour indices with values going from high to low
%
%                   NOTE: x1 and x2 do not need to have the same length
%                   although m and n are usually around 48 - 50
%       
%   xOutline_mc : [96 x n_frames] x values of the contour
%   yOutline_mc : [96 x n_frames] y values of the contour
%   gridSize    : (scalar, normally 50) # of points to use between the minimum and maximum value
%   iFrame      : (scalar) current frame index
%   mn_x_I      : [1 x n_frames] array of minimum x values for all frames
%   mx_x_I      : [1 x n_frames] "       " maximum "           "

X_in_1 = xOutline_mc(x1,iFrame);
X_in_2 = xOutline_mc(x2,iFrame);
Y_in_1 = yOutline_mc(x1,iFrame);
Y_in_2 = yOutline_mc(x2,iFrame);

X_out  = linspace(...
    xOutline_mc(mn_x_I(iFrame),iFrame),...
    xOutline_mc(mx_x_I(iFrame),iFrame),...
    gridSize);

%NOTE: In Matlab the interp1 command has a lot of overhead, so we call the
%real function directly. This requires ordered x values
F = griddedInterpolant(X_in_1,Y_in_1,'linear');
y_interp_1 = F(X_out);

%The values from x2 are going from high to low, so the difference is
%negative. Unfortunately griddedInterpolant requires positve differences,
%so we flip the X and Y vectors
F = griddedInterpolant(X_in_2(end:-1:1),Y_in_2(end:-1:1),'linear');
y_interp_2 = F(X_out);

x_out_all = X_out;

end

function [eccentricity,orientation] = h__getEccentricityAndOrientation(xOutline_mc,yOutline_mc,xRange_all,yRange_all,gridAspectRatio_all,gridSize,eccentricity,orientation,run_mask)

%This is very similar to the way the old code ran. It uses a
%point-in-polygon test for a set of grid points that fill the bounding box
%occupied by the worm. 

%The only functional difference at this point compared to the original code
%is that I've rotated the worms to try and minimize the size of the
%bounding box.

for iFrame = find(run_mask)
    gridAspectRatio = gridAspectRatio_all(iFrame);
    xRange = xRange_all(iFrame);
    yRange = yRange_all(iFrame);
    
    % if the bounding box aspect ratio is skewed, use a smaller number of
    % grid points for the smaller of the two dimensions.  This ensures more
    % uniform sampling of the region inside the polygon.
    if xRange > yRange
        % x size is larger so scale down the number of grid points in the y
        % direction
        wtf1 = linspace(min(xOutline_mc(:,iFrame)), max(xOutline_mc(:,iFrame)), gridSize);
        wtf2 = linspace(min(yOutline_mc(:,iFrame)), max(yOutline_mc(:,iFrame)), round(gridSize / gridAspectRatio));
    else
        % y size is larger so scale down the number of grid points in the x direction
        wtf1 = linspace(min(xOutline_mc(:,iFrame)), max(xOutline_mc(:,iFrame)), round(gridSize * gridAspectRatio));
        wtf2 = linspace(min(yOutline_mc(:,iFrame)), max(yOutline_mc(:,iFrame)), gridSize);
    end
    
    [m,n] = meshgrid( wtf1 , wtf2 );
    
    % get the indices of the points inside of the polygon
    inPointInds = helper__inpolyNew([m(:) n(:)], [xOutline_mc(:,iFrame) yOutline_mc(:,iFrame)]);
    
    %TODO: Send bug report to guy, frame 979 shows that his code does
    %not work all that well ...
    % % %     if iFrame == 979
    % % %        keyboard
    % % %     end
    
    %{
        scatter(xOutline_mc(:,iFrame),yOutline_mc(:,iFrame),'g')
        hold on
        scatter(x,y,'r')
        hold off
        axis equal
    %}
    
    % get the x and y coordinates of the new set of points to be used in calculating eccentricity.
    x = m(inPointInds);
    y = n(inPointInds);
    
    [eccentricity(iFrame),orientation(iFrame)] = h__calculateSingleValues(x,y);
    
end

end

function [eccentricity_s,orientation_s] = h__calculateSingleValues(x,y)

%NOTE: I tried vectorizing this code (all frames at once)
%and it turned out to be slower ...

N = length(x);

% Calculate normalized second central moments for the region.
uxx = sum(x.^2)/N;
uyy = sum(y.^2)/N;
uxy = sum(x.*y)/N;

% Calculate major axis length, minor axis length, and eccentricity.
common          = sqrt((uxx - uyy)^2 + 4*uxy^2);
majorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy + common);
minorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy - common);
eccentricity_s  = 2*sqrt((majorAxisLength/2)^2 - (minorAxisLength/2)^2) / majorAxisLength;

% Calculate orientation.
if (uyy > uxx)
    num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
    den = 2*uxy;
else
    num = 2*uxy;
    den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
end

orientation_s = (180/pi) * atan(num/den);



end

function [cn,on] = helper__inpolyNew(p,node)
%   INPOLY: Point-in-polygon testing.
%
%   Darren Engwirda: 2005-2007
%   Email          : d_engwirda@hotmail.com
%   Last updated   : 23/11/2007 with MATLAB 7.0
%
% Problems or suggestions? Email me.
%
%   JAH: The only change I made was to remove the input checking because
%   the way we call this function will not change.
%
%   Inputs
%   =======================================================================
%   p    [n 2] - points to test if in polygon, each row contains a [x,y]
%                pair
%   node [n 2] - points forming the polygon, NOTE: We assume that the
%                polygon points are linked such that the first and last 
%                points connect

TOL = 1.0e-12;

nnode = size(node,1);
edge  = [(1:nnode-1)' (2:nnode)'; nnode 1]; %This assumes an ordered contour
%input, which is what we provide ...

%% PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n  = size(p,1);
nc = size(edge,1);

% Choose the direction with the biggest range as the "y-coordinate" for the
% test. This should ensure that the sorting is done along the best
% direction for long and skinny problems wrt either the x or y axes.
dxy = max(p,[],1)-min(p,[],1);
if dxy(1)>dxy(2)
    % Flip co-ords if x range is bigger
    p = p(:,[2,1]);
    node = node(:,[2,1]);
end
tol = TOL*min(dxy);

% Sort test points by y-value
[y,i] = sort(p(:,2));
x = p(i,1);

%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cn = false(n,1);     % Because we're dealing with mod(cn,2) we don't have
% to actually increment the crossing number, we can
% just flip a logical at each intersection (faster!)
on = cn;
for k = 1:nc         % Loop through edges
    
    % Nodes in current edge
    n1 = edge(k,1);
    n2 = edge(k,2);
    
    % Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
    %           - also get xmin = min(x1,x2), xmax = max(x1,x2)
    y1 = node(n1,2);
    y2 = node(n2,2);
    if y1<y2
        x1 = node(n1,1);
        x2 = node(n2,1);
    else
        yt = y1;
        y1 = y2;
        y2 = yt;
        x1 = node(n2,1);
        x2 = node(n1,1);
    end
    
    if x1>x2
        xmin = x2;
        xmax = x1;
    else
        xmin = x1;
        xmax = x2;
    end

    % Binary search to find first point with y<=y1 for current edge
    %
    %   SLOW
    %
    if y(1)>=y1
        start = 1;
    elseif y(n)<y1
        start = n+1;
    else
        lower = 1;
        upper = n;
        for j = 1:n
            start = round(0.5*(lower+upper));
            if y(start)<y1
                lower = start;
            elseif y(start-1)<y1
                break;
            else
                upper = start;
            end
        end
    end
    
    % Loop through points
    for j = start:n
        % Check the bounding-box for the edge before doing the intersection
        % test. Take shortcuts wherever possible!
        
        Y = y(j);   % Do the array look-up once & make a temp scalar
        if Y<=y2
            X = x(j);   % Do the array look-up once & make a temp scalar
            if X>=xmin
                if X<=xmax
                    
                    % Check if we're "on" the edge
                    on(j) = on(j) || (abs((y2-Y)*(x1-X)-(y1-Y)*(x2-X))<tol);
                    
                    % Do the actual intersection test
                    if (Y<y2) && ((y2-y1)*(X-x1)<(Y-y1)*(x2-x1))
                        cn(j) = ~cn(j);
                    end
                    
                end
            elseif Y<y2   % Deal with points exactly at vertices
                % Has to cross edge
                cn(j) = ~cn(j);
            end
        else
            % Due to the sorting, no points with >y
            % value need to be checked
            break
        end
    end
    
end

% Re-index to undo the sorting
cn(i) = cn|on;
on(i) = on;

end      % inpoly2()
