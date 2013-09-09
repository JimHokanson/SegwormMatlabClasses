function cleanWorm(obj)
%cleanWorm  Clean up the worm contour by connecting any splits ends.
%
%   cleanWorm(obj)
%
%   Note: the worm's contour is still rough, especially at any split ends.
%         Therefore, index lengths, as opposed to chain-code lengths, are
%         used as the distance metric over the worm's contour.
%
%   Inputs:
%       contour     - the clockwise, circularly-connected worm contour.
%
% wormSegSize -
% The size (in contour points) of a worm segment. Note: The worm's contour
% is roughly divided into 50 segments of musculature (i.e., hinges that
% represent degrees of freedom). Warning: before cleaning, the length of
% the contour can vary significantly: from 1/4 its correct size, if the
% worm is coiled up with its head and tail touching its body, 180 degrees
% apart on the coil; to 2 times its correct size, if the head and tail are
% both split by invaginations that reach 1/4 into its body. Additionally,
% there are various permutations in between these extremes. Therefore, we
% use carefully chosen approximations that are fail-safe to within a large
% margin. Moreover, we use several other tricks to ensure we don't
% incorrectly heal false worm splits (e.g., we check for a sharp concavity
% before joining sharp convexities). But, we remain labile in extreme cases
% (e.g., omega bends where the head and tail are very proximal).
%
%   Output:
%       contour - the cleaned up worm contour.
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

%The following code might be better off if we just did a call to
%computeAngles()
%--------------------------------------------------------------------------
pixels_local  = obj.pixels;
worm_seg_size = size(pixels_local,1)/obj.N_SEGS;

worm_seg_size = round(worm_seg_size);

angles        = seg_worm.cv.circCurvature(pixels_local, worm_seg_size);

% On a small scale, noise causes contour imperfections that shift an angle
% from its correct location. Therefore, blurring angles by averaging them
% with their neighbors can localize them better.
blur_length = ceil(worm_seg_size/2);
mAngles     = seg_worm.util.circConv(angles,[],blur_length);

% Is the worm contour split at the head and/or tail?
% Note: often the head and tail have light colored internals that, when
% binarized, split the head and/or tail into two or more pieces.
% Note 2: We don't use the blurred angles for concavities. Unfortunately,
% blurring can erase high-frequency minima. Moreover, we don't need
% any improvements in localizing these concavities.
%
%   ???? - note, this was a value of 90, but later for head/tail
%   we use 60 for HF and 90 for lf.  This next line is HF, not sure
%   how to reconcile
[~,maxI] = seg_worm.util.peaksCircDist(mAngles, worm_seg_size,true,obj.LF_MAX_ANGLE_LIMIT);
[~,minI] = seg_worm.util.peaksCircDist(mAngles, worm_seg_size,false,obj.MIN_ANGLE_LIMIT);

maxI = maxI(:);
minI = minI(:);
%--------------------------------------------------------------------------
%END OF CODE THAT MIGHT BE REPLACED WITH METHODS CALL
%--------------------------------------------------------------------------





NEAR_SIZE = 2*worm_seg_size; %Used for saying the peaks are too close ...

if length(maxI) > 1
    circ_length_last_to_first_point = (maxI(1) + length(mAngles) - maxI(end));
    maxD = [diff(maxI); circ_length_last_to_first_point];
    
    % Do we have multiple sharp convexities (potential contour splits) that are
    % nearby on the contour and/or, nearby in distance and separated by a sharp
    % concavity?
    %nearScale = .3; % a nearby location on the contour (relative to its size)
    if any(maxD <= NEAR_SIZE) || ~isempty(minI)
        %keyboard
        %Frame 174 of video I was given
        pixels_local = helper__fixInvaginations(obj,maxI,pixels_local,NEAR_SIZE);
    end
end

n_contour_points = size(pixels_local,1);
if n_contour_points > 2
    pixels_local = helper__removeUnecessaryContourPoints(pixels_local);
end

obj.pixels = pixels_local;

end

%
function contour_out = helper__removeUnecessaryContourPoints(contour)
%
%1) For each point, we get the next two points
%2) We remove the following patterns:
%   - after a point, goes back to the previous location
%   - after a point, is on the diogonal from the previous
%
%       see code for details
%
%   Improvements:
%   - This could be sped up by finding problems first then
%   going from problem to problem, instead of point to point

n_contour_points = size(contour,1);

keep_mask = true(1,n_contour_points);
if isequal(contour(1,:), contour(end,:))
    %Remove the first point if it equals the last point
    keep_mask(1) = false;
end

% Remove small overlapping segments and anti alias the contour.
iPoint = 1;
while iPoint <= n_contour_points
    
    % Initialize the next 2 indices.
    %----------------------------------------------------
    if iPoint < n_contour_points - 1
        nextI  = iPoint + 1;
        next2I = iPoint + 2;
    elseif iPoint < size(contour, 1)  % The second index wraps.
        nextI  = size(contour, 1);
        next2I = find(keep_mask,1);
        
        % The are no more kept points.
        if iPoint == next2I
            break
        end
    else  % Both indices wrap.
        
        temp = find(keep_mask,2);
        if length(temp) ~= 2
            break
        end
        
        if any(temp == iPoint)
            break;
        end
        
        nextI  = temp(1);
        next2I = temp(2);
    end
    
    % Remove any overlap.
    %----------------------------------------------------
    dContour = abs(contour(iPoint,:) - contour(next2I,:));
    if all(dContour == 0) %Contour goes back to same location
        %
        %   [ 2 ]
        %   [1,3]  %Remove 1 & 2
        keep_mask(iPoint) = false;
        keep_mask(nextI)  = false;
        
        iPoint = iPoint + 2;
    elseif all(dContour <= 1) % Smooth any stairs.
        %
        %
        %   [2][3]  %<- remove 2
        %   [1]
        keep_mask(nextI) = false;
        
        iPoint = iPoint + 2;
    else %Just advance
        iPoint = iPoint + 1;
    end
end


contour_out = contour(keep_mask,:);
end

%Fix multiple sharp convexities that are nearby in distance and separated
%by a sharp concavity
function contour = helper__fixInvaginations(obj,maxI,contour,NEAR_SIZE)
%
%   JAH NOTE: I didn't look at this code all that closely
%
% Do we have multiple sharp convexities (potential contour splits) that are
% nearby on the contour and/or, nearby in distance and separated by a sharp
% concavity?
%
%
% Connect sharp convexities that are nearby on the contour and/or,
% nearby in distance and separated by a sharp concavity.
% Note: the worm's width is approximately the size of a muscle segment.
% Binarization may yield a split with diagonally-offset, forking
% convexities. Therefore, 2 segments is a good size to bound the
% distance between nearby, split convexities.
% Note 2: the connections are organized as the vector triplet:
% [startContourIndex endContourIndex isWrapping]
% Contour points between startContourIndex and endContourIndex are removed.

%??? - does this mean invaginations exist if this code is called?
%Where is there a yes or no

conns  = zeros(length(maxI), 4); % the connections (pre-allocate memory)
connsI = 1; % the current index for connections
for i = 1:(length(maxI) - 1);
    
    % Are there any sharp convexities nearby?
    for j = (i + 1):length(maxI)
        if sqrt(sum((contour(maxI(i),:) - contour(maxI(j),:)) .^ 2)) <= NEAR_SIZE
            
            % Which side is shorter?
            % Side1 is continuous and goes from start (iI) to end (jI)
            % in positive, index increments.
            % Side2 wraps and always goes from start (iI) to end (jI)
            % in negative, index increments.
            cLength = size(contour, 1);
            iI = maxI(i);
            jI = maxI(j);
            dSide1 = jI - iI;
            dSide2 = iI + cLength - jI;
            
            % The continuous side is shorter.
            if dSide1 < dSide2 % && dSide1 / cLength < nearScale
                
                % Is the convexity nearby on the contour.
                if dSide1 <= NEAR_SIZE
                    conns(connsI,:) = [iI jI 0 dSide1];
                    connsI = connsI + 1;
                    
                    % Is there a concavity separating us on our shorter,
                    % continuous side?
                else
                    for k = 1:length(minI)
                        if minI(k) > iI && minI(k) < jI
                            conns(connsI,:) = [iI jI 0 dSide1];
                            connsI = connsI + 1;
                            break;
                        end
                    end
                end
                
                % The wrapping side is shorter so check it instead.
            else %if dSide2 / cLength < nearScale
                
                % Is the convexity nearby on the contour.
                if dSide2 <= NEAR_SIZE
                    conns(connsI,:) = [jI iI 1 dSide2];
                    connsI = connsI + 1;
                    
                    % Is there a concavity separating us on our shorter,
                    % wrapping side?
                else
                    for k = 1:length(minI)
                        if minI(k) < iI || minI(k) > jI
                            conns(connsI,:) = [jI iI 1 dSide2];
                            connsI = connsI + 1;
                            break;
                        end
                    end
                end
                % ['Connecting ' num2str(maxI(i)) ' and ' num2str(maxI(j)) ...
                %     ' with ' num2str(minI(k)) ' in between.']
            end
        end
    end
end

% Collapse any extra memory.
conns(connsI:end,:) = [];

% Clean up the contour.
if ~isempty(conns)
    
    % Sort the connections.
    conns = sortrows(conns, 4);
    
    % Connect the peaks until there are at least 2 left.
    peaks = zeros(length(maxI),1);
    numPeaks = length(maxI);
    if numPeaks > 2
        peaks(1:2,1) = conns(1,1:2); % connect the peaks
        peaks(1:2,2) = 1; % label the new, unique peak connection
        j = 3; % the peaks index
        label = 2; % the unique peak label index
        numPeaks = numPeaks - 1; % the number of unique peaks
        i = 2; % the conns index
        while numPeaks > 2 && i <= size(conns, 1)
            
            % Are either of the peaks new?
            peak1 = peaks(peaks(1:(j - 1),1) == conns(i,1),:);
            peak2 = peaks(peaks(1:(j - 1),1) == conns(i,2),:);
            
            % Both peaks are new.
            if isempty(peak1)
                if isempty(peak2)
                    peaks(j:(j + 1),1) = conns(i,1:2);
                    peaks(j:(j + 1),2) = label;
                    j = j + 2;
                    label = label + 1;
                    
                    % The first peak is new.
                else
                    peaks(j,:) = [conns(i,1) peak2(2)];
                    j = j + 1;
                end
                
                % We lost a peak to the connection.
                numPeaks = numPeaks - 1;
                
                % The second peak is new.
            elseif isempty(peak2)
                peaks(j,:) = [conns(i,2) peak1(2)];
                j = j + 1;
                
                % We lost a peak to the connection.
                numPeaks = numPeaks - 1;
                
                % Relabel the second peak and its connections.
            elseif peak1(2) < peak2(2)
                peaks(peaks(1:(j - 1),2) == peak2(2),2) = peak1(2);
                
                % We lost a peak to the connection.
                numPeaks = numPeaks - 1;
                
                % Relabel the first peak and its connections.
            elseif peak1(2) > peak2(2)
                peaks(peaks(1:(j - 1),2) == peak1(2),2) = peak2(2);
                
                % We lost a peak to the connection.
                numPeaks = numPeaks - 1;
            end
            
            % Advance.
            i = i + 1;
        end
        conns(i:end,:) = [];
    end
    
    % Connect the connections.
    prevConnsSize = size(conns, 1);
    newConnsI = 1; % the current index for new connections
    while newConnsI < prevConnsSize
        newConns = zeros(size(conns, 1), 3); % the new connections (pre-allocate memory)
        newConnsI = 1;
        for i = 1:size(conns, 1)
            connected = false; % have we made any connections?
            for j = (i + 1):size(conns, 1)
                
                % Are both connections continuous?
                if ~conns(i,3)
                    if ~conns(j,3)
                        
                        % Does connection j intersect i?
                        if conns(i,2) - conns(i,1) >= conns(j,2) - conns(j,1)
                            if (conns(i,1) <= conns(j,1) && conns(i,2) >= conns(j,1)) ...
                                    || (conns(i,1) <= conns(j,2) && conns(i,2) >= conns(j,2))
                                
                                % Take the union of connections i and j.
                                newConns(newConnsI,:) = ...
                                    [ min(conns(i,1), conns(j,1)) ...
                                    max(conns(i,2), conns(j,2)) ...
                                    0 ];
                                newConnsI = newConnsI + 1;
                                connected = true;
                            end
                            
                            % Does connection i intersect j?
                        else
                            if (conns(i,1) >= conns(j,1) && conns(i,1) <= conns(j,2)) ...
                                    || (conns(i,2) >= conns(j,1) && conns(i,2) <= conns(j,2))
                                
                                % Take the union of connections i and j.
                                newConns(newConnsI,:) = ...
                                    [ min(conns(i,1), conns(j,1)) ...
                                    max(conns(i,2), conns(j,2)) ...
                                    0 ];
                                newConnsI = newConnsI + 1;
                                connected = true;
                            end
                        end
                        
                        % Connection j wraps.
                    else
                        
                        % Add connection i to the beginning of j.
                        justConnected = false; % did we just connect?
                        if conns(i,2) >= conns(j,1)
                            newConns(newConnsI,:) = ...
                                [ min(conns(i,1), conns(j,1))
                                conns(j,2)
                                1 ];
                            newConnsI = newConnsI + 1;
                            connected = true;
                            justConnected = true;
                        end
                        
                        % Add connection i to the end of j.
                        if conns(i,1) <= conns(j,2)
                            if justConnected
                                newConns(newConnsI - 1,2) = ...
                                    max(conns(i,2), conns(j,2));
                            else
                                newConns(newConnsI,:) = ...
                                    [ conns(j,1)
                                    max(conns(i,2), conns(j,2))
                                    1 ];
                                newConnsI = newConnsI + 1;
                                connected = true;
                            end
                        end
                    end
                    
                    % Are both connections wrapping?
                else
                    if conns(j,3)
                        
                        % Take the union of connections i and j.
                        newConns(newConnsI,:) = ...
                            [ min(conns(i,1), conns(j,1)) ...
                            max(conns(i,2), conns(j,2)) ...
                            1 ];
                        newConnsI = newConnsI + 1;
                        connected = true;
                        
                        % Connection j is continuous.
                    else
                        
                        % Add connection j to the beginning of i.
                        justConnected = false; % did we just connect?
                        if conns(i,1) <= conns(j,2)
                            newConns(newConnsI,:) = ...
                                [ min(conns(i,1), conns(j,1))
                                conns(i,2)
                                1 ];
                            newConnsI = newConnsI + 1;
                            connected = true;
                            justConnected = true;
                        end
                        
                        % Add connection j to the end of i.
                        if conns(i,2) >= conns(j,1)
                            if justConnected
                                newConns(newConnsI - 1,2) = ...
                                    max(conns(i,2), conns(j,2));
                            else
                                newConns(newConnsI,:) = ...
                                    [ conns(i,1)
                                    max(conns(i,2), conns(j,2))
                                    1 ];
                                newConnsI = newConnsI + 1;
                                connected = true;
                            end
                        end
                    end
                end
            end
            
            % Add the connection.
            if ~connected
                newConns(newConnsI,:) = conns(i,1:3);
                newConnsI = newConnsI + 1;
            end
        end
        
        % Collapse any extra memory.
        newConns(newConnsI:end,:) = [];
        
        % Have we made any new connections?
        prevConnsSize = size(conns, 1);
        conns = newConns;
    end
    
    % Connect the contour splits.
    for i = 1:size(conns, 1)
        
        % Connect the continuous contour split.
        if ~conns(i,3)
            minI = conns(i,1);
            maxI = conns(i,2);
            minP = contour(minI,:);
            maxP = contour(maxI,:);
            points = maxI - minI + 1;
            contour(minI:maxI,1) = round(linspace(minP(1), maxP(1), points));
            contour(minI:maxI,2) = round(linspace(minP(2), maxP(2), points));
            
            % Connect the wrapping contour split.
        else
            minI   = conns(i,2);
            maxI   = conns(i,1);
            minP   = contour(minI,:);
            maxP   = contour(maxI,:);
            points = minI + size(contour, 1) - maxI + 1;
            interPoints = [];
            interPoints(:,1) = linspace(maxP(1), minP(1), points);
            interPoints(:,2) = linspace(maxP(2), minP(2), points);
            contour(maxI:end,:) = round(interPoints(1:(end - minI),:));
            contour(1:minI,:)  = round(interPoints((end - minI + 1):end,:));
        end
    end
    
    % Clean up the contour.
    %temp_pixels = obj.pixels;
    obj.pixels = obj.cleanContour(obj.pixels);
    %final_pixels = obj.pixels;
    
%     cplot = seg_worm.video.pixel_list_image;
%     cplot.addList('original',[0 0 1],temp_pixels)
%     cplot.addList('final',[0 1 0],final_pixels)
%     cplot.plot();
%     keyboard
    
    %contour = obj.pixels();
    %contour = cleanContour(contour);
end




end