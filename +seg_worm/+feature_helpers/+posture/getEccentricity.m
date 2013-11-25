function [eccentricity, orientation] = getEccentricity(xOutline, yOutline, gridSize)
% GETPOLYECCENTRICITY  Given x and y coordinates of the outline of a region
% of interest, fill the outline with a grid of evenly spaced points and use
% these to calculate the eccentricity and orientation of the equivalent
% ellipse.
%
%
%   I think this could be sped up a bit ...
%
% Input:
%   xOutline - The x coordinates of the outline of the region of interest
%   yOutline - The y coordinates of the outline of the region of interest
%   gridSize - The side length of grid to be used to fill in the region of
%              interest.  50 gives reasonable values quickly for our data.
%
% Output:
%   eccentricity - The eccentricity of the equivalent ellipse of the region
%                  of interest.
%   orientation  - The orientation angle of the equivalent ellipse of the
%                  region of interest.
%
% NOTE: this function uses the eccentricity code from regionprops
% essentially unchanged.
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

MIDDLE_INDEX = 49;



%OLD CODE
%--------------------------------------------------------------------------

%{
xOutline_mc = bsxfun(@minus,xOutline,mean(xOutline,1)); %_mc - mean centered
yOutline_mc = bsxfun(@minus,yOutline,mean(yOutline,1));

xRange_all = max(xOutline_mc,[],1) - min(xOutline_mc,[],1);
yRange_all = max(yOutline_mc,[],1) - min(yOutline_mc,[],1);
gridAspectRatio_all = xRange_all./yRange_all;

tic
[eccentricity1,orientation1] = h__getEccentricityAndOrientationOld(xOutline_mc,yOutline_mc,xRange_all,yRange_all,gridAspectRatio_all,gridSize);
toc
%}

%Rough time breakdown
%--------------------------------------
%1) Setup - 0.38
%2) y value calculations: 0.004
%3) simple worm counting: 0.555
%4) remaining grid points: 0.440
%
%   3035 - simple
%   325  - uses grid to calculate
%

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
%   points on the grid outside of the worm).
%
%   2) This allows us to hardcode only looking for "simple" worms 
%   (see description below) that vary in the x-direction.
%

y = yOutline_mc(MIDDLE_INDEX,:) - yOutline_mc(1,:);
x = xOutline_mc(MIDDLE_INDEX,:) - xOutline_mc(1,:);

rot = atan2(y,x);

%I expanded the rotation matrix to allow processing all frames at once
%
%   i.e. rather than R*X (Matrix multiplication)
%
%   I have r(1)*X(1) + r(2)*X(2) etc, but where X(1), X(2), etc is really
%   a vector of values, not just a singular value like you would need for
%   R*X
%
xOutline_rot = bsxfun(@times,xOutline_mc,cos(rot)) + bsxfun(@times,yOutline_mc,sin(rot));
yOutline_rot = bsxfun(@times,xOutline_mc,-sin(rot)) + bsxfun(@times,yOutline_mc,cos(rot));

%TODO: replace refs of _mc with _rot, this is a temporary fix
xOutline_mc = xOutline_rot; 
yOutline_mc = yOutline_rot;

%{

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
%   Simple worms have two sides that go from - to +, with no bending
%   backwards. If we find a simple worm then we can interpolate the y
%   values at the x-grid locations. This removes the need to sort the x
%   values in the grid with respect to the x-values of the contour. Once
%   we have the x-values aligned, we round the y-values to grid points.
%   Doing this allows us to simply count the points off that are between
%   the two y-values.
%
%   For the y-values, consider the points 13 and 33, and grid values
%   at every 5, starting at 11. We take 13 and round up to the nearest grid
%   point, which would be 16. We take 33 and round down to the nearest grid
%   point, which would be floor((33 - 11)/5) => 20, 20 + 11 => 31. So now
%   that we have 16 and 31 the points between 13 and 33 are:
%   16:5:31 or 16, 21, 26, 31
%
[mx_x,mx_x_I] = max(xOutline_rot,[],1);
[mn_x,mn_x_I] = min(xOutline_rot,[],1);

min_first_mask = mn_x_I < mx_x_I;
min_last_mask  = mx_x_I < mn_x_I; %NOTE: We are excluding equals ...


mx_y = max(yOutline_mc,[],1);
mn_y = min(yOutline_mc,[],1);

d = [diff(xOutline_mc,1,1); xOutline_mc(1,:) - xOutline_mc(end,:)];

n_frames = size(xOutline,2);
n_c = size(xOutline,1);

y_interp_1 = NaN(gridSize,n_frames);
y_interp_2 = NaN(gridSize,n_frames);
x_out_all  = NaN(gridSize,n_frames);


%TODO: Could remove linspace

%--------------------------------------------------------------------------
%x1: from negative to postive x
%x2: form positive to negative x - needs to be reversed during interpolation
for iFrame = find(min_first_mask)
    
    x1 = mn_x_I(iFrame):(mx_x_I(iFrame));
    
    if all(d(x1(1:end-1),iFrame) > 0)
        
        x2 = [mx_x_I(iFrame):n_c 1:(mn_x_I(iFrame))];
        
        if all(d(x2(1:end-1),iFrame) < 0)
            [y_interp_1(:,iFrame),y_interp_2(:,iFrame),x_out_all(:,iFrame)] = ...
                h__getInterpValues(x1,x2,xOutline_mc,yOutline_mc,gridSize,iFrame,mn_x_I,mx_x_I);
        end
    end
end

%--------------------------------------------------------------------------
for iFrame = find(min_last_mask)
    
    x2 = mx_x_I(iFrame):(mn_x_I(iFrame));
    
    if all(d(x2(1:end-1),iFrame) < 0)
        
        x1 = [mn_x_I(iFrame):n_c 1:(mx_x_I(iFrame))];
        
        if all(d(x1(1:end-1),iFrame) > 0)
            [y_interp_1(:,iFrame),y_interp_2(:,iFrame),x_out_all(:,iFrame)] = ...
                h__getInterpValues(x1,x2,xOutline_mc,yOutline_mc,gridSize,iFrame,mn_x_I,mx_x_I);
        end
    end
end

%--------------------------------------------------------------------------

%TODO: Now that we have have points at exactly the same X, we need
%to determine the y values, we might be able to do this by rounding

TOL = 1.0e-12;
%TODO: 

% max_y_interp = max(y_interp_2,[],1);
min_y_interp = min(y_interp_1,[],1);


% nY_all = round(gridSize./gridAspectRatio_all);
% d_unit = (max_y_interp - min_y_interp)./(nY_all - 1);

%NOTE: d_unit should be the same as the x_unit, otherwise we are biased
d_unit = x_out_all(2,:) - x_out_all(1,:);


x_all = zeros(1,gridSize*gridSize);
y_all = zeros(1,gridSize*gridSize);

eccentricity = NaN(1,n_frames);
orientation  = NaN(1,n_frames);

simple_worm_mask = any(x_out_all,1);

d_unit = d_unit(simple_worm_mask);
min_y_interp = min_y_interp(simple_worm_mask);

%This step takes a long time ...
%Can we speed it up????
%
%Could vectorize rounding

%TODO: I need to ensure that y2 > y1, if not, flip
%NOTE: beginning and end will be a tossup due to rounding error, any 
%point in the middle should be fine for testing

y_1_rounded = bsxfun(@minus,y_interp_1(:,simple_worm_mask),min_y_interp)-TOL;
y_2_rounded = bsxfun(@minus,y_interp_2(:,simple_worm_mask),min_y_interp)-TOL;
y_1_rounded = bsxfun(@times,ceil(bsxfun(@rdivide,y_1_rounded,d_unit)),d_unit);
y_2_rounded = bsxfun(@times,floor(bsxfun(@rdivide,y_2_rounded,d_unit)),d_unit);
y_1_rounded = bsxfun(@plus,y_1_rounded,min_y_interp);
y_2_rounded = bsxfun(@plus,y_2_rounded,min_y_interp);

cur_simple = 0;
for iFrame = find(simple_worm_mask)
    count = 0;
    cur_simple = cur_simple + 1;
    
    cur_d_unit = d_unit(cur_simple);
    
    for iIndex = 1:gridSize
        temp = y_1_rounded(iIndex,cur_simple):cur_d_unit:y_2_rounded(iIndex,cur_simple);
        y_all(count+1:count+length(temp)) = temp;
        x_all(count+1:count+length(temp)) = x_out_all(iIndex,iFrame);
        count = count + length(temp);
    end

    [eccentricity(iFrame),orientation(iFrame)] = h__calculateSingleValues(x_all(1:count),y_all(1:count));    
end

%{
orientation2 = orientation;
eccentricity2 = eccentricity;
%}

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
%TODO: Pass in min and max, and compute these values inside
%the function 
xRange_all = mx_x - mn_x;
yRange_all = mx_y - mn_y;
gridAspectRatio_all = xRange_all./yRange_all;

run_mask = ~simple_worm_mask & ~isnan(gridAspectRatio_all);

[eccentricity,orientation] = h__getEccentricityAndOrientation(xOutline_mc,yOutline_mc,xRange_all,yRange_all,gridAspectRatio_all,gridSize,eccentricity,orientation,run_mask);

%Fix the orientation ...
%--------------------------------------------------------------------------
orientation3 = orientation + rot*180/pi;
orientation3(orientation3 > 90) = orientation3(orientation3 > 90) - 180;
orientation3(orientation3 < -90) = orientation3(orientation3 < -90) + 180;

orientation = orientation3;

%keyboard


%{
subplot(2,1,1)
plot(eccentricity - eccentricity1)
subplot(2,1,2)
plot(orientation - orientation1)

%}

end

function [y_interp_1,y_interp_2,x_out_all] = h__getInterpValues(x1,x2,xOutline_mc,yOutline_mc,gridSize,iFrame,mn_x_I,mx_x_I)

%Does interpolation for the two halves
%TODO: Define inputs ...

X_in_1 = xOutline_mc(x1,iFrame);
X_in_2 = xOutline_mc(x2,iFrame);
Y_in_1 = yOutline_mc(x1,iFrame);
Y_in_2 = yOutline_mc(x2,iFrame);
X_out  = linspace(...
    xOutline_mc(mn_x_I(iFrame),iFrame),...
    xOutline_mc(mx_x_I(iFrame),iFrame),...
    gridSize);

F = griddedInterpolant(X_in_1,Y_in_1,'linear');
y_interp_1 = F(X_out);

F = griddedInterpolant(X_in_2(end:-1:1),Y_in_2(end:-1:1),'linear');
y_interp_2 = F(X_out);

x_out_all = X_out;

end

function [eccentricity,orientation] = h__getEccentricityAndOrientation(xOutline_mc,yOutline_mc,xRange_all,yRange_all,gridAspectRatio_all,gridSize,eccentricity,orientation,run_mask)

%The goal of this code is to go from a contour to 'pixels' that are inside
%the contour. Using these pixels as sample points we can fit an ellipse to
%the data ...
%==========================================================================

for iFrame = find(run_mask)
    gridAspectRatio = gridAspectRatio_all(iFrame);
    xRange = xRange_all(iFrame);
    yRange = yRange_all(iFrame);
    %NOTE: The data needs to be sampled evenly, otherwise it will be biased ...
    
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
    %Are the points
    
    %TODO: Send bug report to guy
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


function [eccentricity,orientation] = h__getEccentricityAndOrientationOld(xOutline_mc,yOutline_mc,xRange_all,yRange_all,gridAspectRatio_all,gridSize)

%The goal of this code is to go from a contour to 'pixels' that are inside
%the contour. Using these pixels as sample points we can fit an ellipse to
%the data ...
%==========================================================================


n_frames = size(xOutline_mc,2);
eccentricity = NaN(1,n_frames);
orientation  = NaN(1,n_frames);

I_use = find(~isnan(gridAspectRatio_all));
I_use = I_use(:)'; %We need this to be a row vector for iterating over ...

for iFrame = I_use
    gridAspectRatio = gridAspectRatio_all(iFrame);
    xRange = xRange_all(iFrame);
    yRange = yRange_all(iFrame);
    %NOTE: The data needs to be sampled evenly, otherwise it will be biased ...
    
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
    inPointInds = helper__inpoly([m(:) n(:)], [xOutline_mc(:,iFrame) yOutline_mc(:,iFrame)]);
    %Are the points
    
    
    % get the x and y coordinates of the new set of points to be used in calculating eccentricity.
    x = m(inPointInds);
    y = n(inPointInds);
    %This basically goes back to
    
    %==========================================================================
    
    
    N = length(x);
    
    % Calculate normalized second central moments for the region.
    uxx = sum(x.^2)/N;
    uyy = sum(y.^2)/N;
    uxy = sum(x.*y)/N;
    
    % Calculate major axis length, minor axis length, and eccentricity.
    common          = sqrt((uxx - uyy)^2 + 4*uxy^2);
    majorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy + common);
    minorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy - common);
    eccentricity(iFrame) = 2*sqrt((majorAxisLength/2)^2 - (minorAxisLength/2)^2) / majorAxisLength;
    
    % Calculate orientation.
    if (uyy > uxx)
        num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
        den = 2*uxy;
    else
        num = 2*uxy;
        den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
    end
    
    if (num == 0) && (den == 0)
        orientation(iFrame) = 0;
    else
        orientation(iFrame) = (180/pi) * atan(num/den);
    end
end

end

function [cn,on] = helper__inpoly(p,node,edge,TOL)
%  INPOLY: Point-in-polygon testing.
%
% Determine whether a series of points lie within the bounds of a polygon
% in the 2D plane. General non-convex, multiply-connected polygonal
% regions can be handled.
%
% SHORT SYNTAX:
%
%   in = inpoly(p,node);
%
%   p   : The points to be tested as an Nx2 array [x1 y1; x2 y2; etc].
%   node: The vertices of the polygon as an Mx2 array [X1 Y1; X2 Y2; etc].
%         The standard syntax assumes that the vertices are specified in
%         consecutive order.
%
%   in  : An Nx1 logical array with IN(i) = TRUE if P(i,:) lies within the
%         region.
%
% LONG SYNTAX:
%
%  [in,on] = inpoly(p,node,edge);
%
%  edge: An Mx2 array of polygon edges, specified as connections between
%        the vertices in NODE: [n1 n2; n3 n4; etc]. The vertices in NODE
%        do not need to be specified in connsecutive order when using the
%        extended syntax.
%
%  on  : An Nx1 logical array with ON(i) = TRUE if P(i,:) lies on a
%        polygon edge. (A tolerance is used to deal with numerical
%        precision, so that points within a distance of
%        eps^0.8*norm(node(:),inf) from a polygon edge are considered "on"
%        the edge.
%
% EXAMPLE:
%
%   polydemo;       % Will run a few examples
%
% See also INPOLYGON

% The algorithm is based on the crossing number test, which counts the
% number of times a line that extends from each point past the right-most
% region of the polygon intersects with a polygon edge. Points with odd
% counts are inside. A simple implementation of this method requires each
% wall intersection be checked for each point, resulting in an O(N*M)
% operation count.
%
% This implementation does better in 2 ways:
%
%   1. The test points are sorted by y-value and a binary search is used to
%      find the first point in the list that has a chance of intersecting
%      with a given wall. The sorted list is also used to determine when we
%      have reached the last point in the list that has a chance of
%      intersection. This means that in general only a small portion of
%      points are checked for each wall, rather than the whole set.
%
%   2. The intersection test is simplified by first checking against the
%      bounding box for a given wall segment. Checking against the bbox is
%      an inexpensive alternative to the full intersection test and allows
%      us to take a number of shortcuts, minimising the number of times the
%      full test needs to be done.
%
%   Darren Engwirda: 2005-2007
%   Email          : d_engwirda@hotmail.com
%   Last updated   : 23/11/2007 with MATLAB 7.0
%
% Problems or suggestions? Email me.

%% ERROR CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
    TOL = 1.0e-12;
    if nargin<3
        edge = [];
        if nargin<2
            error('Insufficient inputs');
        end
    end
end

nnode = size(node,1);
if isempty(edge)                                                           % Build edge if not passed
    edge = [(1:nnode-1)' (2:nnode)'; nnode 1];
end
if size(p,2)~=2
    error('P must be an Nx2 array.');
end
if size(node,2)~=2
    error('NODE must be an Mx2 array.');
end
if size(edge,2)~=2
    error('EDGE must be an Mx2 array.');
end
if max(edge(:))>nnode || any(edge(:)<1)
    error('Invalid EDGE.');
end

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

end      % inpoly()

function [cn,on] = helper__inpolyNew(p,node)
%  INPOLY: Point-in-polygon testing.
%
%   Darren Engwirda: 2005-2007
%   Email          : d_engwirda@hotmail.com
%   Last updated   : 23/11/2007 with MATLAB 7.0
%
% Problems or suggestions? Email me.

%p    [n 2] - points to test if in polygon
%node [n 2] - points forming the polygon

TOL = 1.0e-12;

nnode = size(node,1);
edge = [(1:nnode-1)' (2:nnode)'; nnode 1]; %This assumes an ordered contour
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

% % % % is_ordered_mask = node(edge(:,1),2) < node(edge(:,2),2);
% % % %
% % % % y1a = zeros(1,nc);
% % % % y2a = zeros(1,nc);
% % % % x1a = zeros(1,nc);
% % % % x2a = zeros(1,nc);
% % % %
% % % % y1a(is_ordered_mask)  = node(edge(is_ordered_mask,1),2);
% % % % y1a(~is_ordered_mask) = node(edge(~is_ordered_mask,2),2);
% % % % y2a(is_ordered_mask)  = node(edge(is_ordered_mask,2),2);
% % % % y2a(~is_ordered_mask) = node(edge(~is_ordered_mask,1),2);
% % % %
% % % % x1a(is_ordered_mask)  = node(edge(is_ordered_mask,1),1);
% % % % x1a(~is_ordered_mask) = node(edge(~is_ordered_mask,2),1);
% % % % x2a(is_ordered_mask)  = node(edge(is_ordered_mask,2),1);
% % % % x2a(~is_ordered_mask) = node(edge(~is_ordered_mask,1),1);
% % % %
% % % % xmina = min(x1a,x2a);
% % % % xmaxa = max(x1a,x2a);
% % % %
% % % % %EDGES(k) <= X(i) < EDGES(k+1)
% % % % %N = histc(X,EDGES)
% % % %
% % % % [~,bin1] = histc(y1a,y);
% % % % [~,bin2] = histc(y2a,y);


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
    
    % % % %     x1 = x1a(k);
    % % % %     x2 = x2a(k);
    % % % %     y1 = y1a(k);
    % % % %     y2 = y2a(k);
    % % % %     xmin = xmina(k);
    % % % %     xmax = xmaxa(k);
    
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
