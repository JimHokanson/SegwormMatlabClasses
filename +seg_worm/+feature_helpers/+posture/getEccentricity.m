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

%center the worm at its centroid
xOutline = bsxfun(@minus,xOutline,mean(xOutline,1));
yOutline = bsxfun(@minus,yOutline,mean(yOutline,1));

% make a grid of points that covers the worm area based on the max and min values of xOutline and yOutline.
xRange_all = max(xOutline,[],1) - min(xOutline,[],1);
yRange_all = max(yOutline,[],1) - min(yOutline,[],1);
gridAspectRatio_all = xRange_all./yRange_all;


%The code below will need to change for multiple frames at once ...
%--------------------------------------------------------------------------

%The goal of this code is to go from a contour to 'pixels' that are inside
%the contour. Using these pixels as sample points we can fit an ellipse to
%the data ...
%==========================================================================

n_frames = size(xOutline,2);
eccentricity = NaN(1,n_frames);
orientation  = NaN(1,n_frames);

I_use = find(~isnan(gridAspectRatio_all));
I_use = I_use(:)';


%OLD VERSION OF CODE


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
        wtf1 = linspace(min(xOutline(:,iFrame)), max(xOutline(:,iFrame)), gridSize);
        wtf2 = linspace(min(yOutline(:,iFrame)), max(yOutline(:,iFrame)), round(gridSize / gridAspectRatio));
    else
        % y size is larger so scale down the number of grid points in the x direction
        wtf1 = linspace(min(xOutline(:,iFrame)), max(xOutline(:,iFrame)), round(gridSize * gridAspectRatio));
        wtf2 = linspace(min(yOutline(:,iFrame)), max(yOutline(:,iFrame)), gridSize);
    end
    
    [m,n] = meshgrid( wtf1 , wtf2 );
    
    % get the indices of the points inside of the polygon
    inPointInds = helper__inpoly([m(:) n(:)], [xOutline(:,iFrame) yOutline(:,iFrame)]);
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
    eccentricity(iFrame)    = 2*sqrt((majorAxisLength/2)^2 - (minorAxisLength/2)^2) / majorAxisLength;
    
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




%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================


%{
tic

orientation2  = NaN(1,n_frames);
uxx2 = NaN(1,n_frames);
uyy2 = NaN(1,n_frames);
uxy2 = NaN(1,n_frames);

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
    wtf1 = linspace(min(xOutline(:,iFrame)), max(xOutline(:,iFrame)), gridSize);
    wtf2 = linspace(min(yOutline(:,iFrame)), max(yOutline(:,iFrame)), round(gridSize / gridAspectRatio));
else
    % y size is larger so scale down the number of grid points in the x direction
    wtf1 = linspace(min(xOutline(:,iFrame)), max(xOutline(:,iFrame)), round(gridSize * gridAspectRatio));
    wtf2 = linspace(min(yOutline(:,iFrame)), max(yOutline(:,iFrame)), gridSize);
end

[m,n] = meshgrid( wtf1 , wtf2 );

% get the indices of the points inside of the polygon
inPointInds = helper__inpoly2([m(:) n(:)], [xOutline(:,iFrame) yOutline(:,iFrame)]);

%inPointInds = inpolygon(m(:),n(:),xOutline(:,iFrame),yOutline(:,iFrame));

%Are the points


% get the x and y coordinates of the new set of points to be used in calculating eccentricity.
x = m(inPointInds);
y = n(inPointInds);
%This basically goes back to

%==========================================================================


N = length(x);

% Calculate normalized second central moments for the region.
uxx2(iFrame) = sum(x.^2)/N;
uyy2(iFrame)  = sum(y.^2)/N;
uxy2(iFrame)  = sum(x.*y)/N;
end

common          = sqrt((uxx2 - uyy2).^2 + 4*uxy2.^2);
majorAxisLength = 2*sqrt(2).*sqrt(uxx2 + uyy2 + common);
minorAxisLength = 2*sqrt(2).*sqrt(uxx2 + uyy2 - common);
eccentricity2   = 2*sqrt((majorAxisLength./2).^2 - (minorAxisLength./2).^2) ./ majorAxisLength;


toc

% % % % Calculate major axis length, minor axis length, and eccentricity.
% % % common          = sqrt((uxx - uyy)^2 + 4*uxy^2);
% % % majorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy + common);
% % % minorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy - common);
% % % eccentricity(iFrame)    = 2*sqrt((majorAxisLength/2)^2 - (minorAxisLength./2).^2) / majorAxisLength;
% % %
% % % % Calculate orientation.
% % % if (uyy > uxx)
% % %     num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
% % %     den = 2*uxy;
% % % else
% % %     num = 2*uxy;
% % %     den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
% % % end
% % %
% % % if (num == 0) && (den == 0)
% % %     orientation(iFrame) = 0;
% % % else
% % %     orientation(iFrame) = (180/pi) * atan(num/den);
% % % end
% % % end
% % % toc

%}


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

function [cn,on] = helper__inpoly2(p,node)
%  INPOLY: Point-in-polygon testing.

%p [n 2]
%node [n 2]

merge = [p; node];


[~,Ix] = sort(merge(:,1));

I1 = find(Ix > size(p,1));

node_indices = zeros(1,size(node,1));

node_indices(Ix(I1)-size(p,1)) = I1;

%For each node, know it's location in Ix, Iy


%% ERROR CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TOL   = 1.0e-12;
nnode = size(node,1);
edge  = [(1:nnode-1)' (2:nnode)'; nnode 1];

edge_indices = node_indices(edge);




%edge, start node, end node

%% PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n  = size(p,1);
nc = size(edge,1);

% Choose the direction with the biggest range as the "y-coordinate" for the
% test. This should ensure that the sorting is done along the best
% direction for long and skinny problems wrt either the x or y axes.
dxy = max(p,[],1) - min(p,[],1);
if dxy(1) > dxy(2)
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

%cn - when point crosses an edge, mark true
%when it crosses another edge, mark false (left the polygon)

on = cn;
for k = 1:nc         % Loop through edges
    
    % Nodes in current edge
    n1 = edge(k,1);
    n2 = edge(k,2);
    
    % Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
    %           - also get xmin = min(x1,x2), xmax = max(x1,x2)
    y1 = node(n1,2);
    y2 = node(n2,2);
    if y1 < y2
        x1 = node(n1,1);
        x2 = node(n2,1);
    else
        yt = y1;
        y1 = y2;
        y2 = yt;
        x1 = node(n2,1);
        x2 = node(n1,1);
    end
    
    if x1 > x2
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
        if Y <= y2
            X = x(j);   % Do the array look-up once & make a temp scalar
            if X >= xmin
                if X <= xmax
                    
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
