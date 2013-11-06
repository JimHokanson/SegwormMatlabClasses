function indices = chainCodeLength2Index(lengths, cc_lengths)
%chainCodeLength2Index  Translate a length into an index. The index
%   represents the numerically-closest element to the desired length in
%   an ascending array of chain code lengths.
%
%   indices = seg_worm.cv.chainCodeLength2Index(lengths, cc_lengths)
%
%   Inputs:
%       lengths    - the lengths to translate into indices
%       cc_lengths - an ascending array of chain code lengths
%                          Note: the chain code lengths must increase at
%                          every successive index
%
%
%   ??? Is this for the cc_lengths or the skeleton cc_lengths?
%
%   The code currently assumes circular chain code lengths
%
%   Output:
%       indices - the indices for the elements closest to the desired lengths
%
%   See also:
%   COMPUTECHAINCODELENGTHS, 
%   CIRCCOMPUTECHAINCODELENGTHS
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

persistent cc_input p_F

% Check the lengths.
% Note: circular chain-code lengths are minimally bounded at 0 and
% maximally bounded at the first + last lengths.
if any(lengths < 0)
    error('chainCodeLength2Index:TooShort','The lengths cannot be negative');
elseif any(lengths > cc_lengths(end))
    error('chainCodeLength2Index:TooLong',...
        'The lengths cannot be greater than %g',cc_lengths(end));
end

if isequal(cc_input,cc_lengths)
    F = p_F;
else
    F = griddedInterpolant(cc_lengths,1:length(cc_lengths),'linear');
    p_F = F;
    cc_input = cc_lengths;
end

indices = round(F(lengths));

%This assumes a circular chain code length ...
indices(indices == 0) = length(cc_lengths);
indices(indices == length(cc_lengths)+1) = 1;


end

function indices = helper__oldCode(lengths, cc_lengths)

% Go through the lengths.
indices = zeros(size(lengths));
for i = 1:numel(lengths)

    % Is the length too small?
    if lengths(i) < cc_lengths(1)
        
        % Find the closest index.
        if lengths(i) / cc_lengths(1) < .5
            indices(i) = length(cc_lengths);
        else
            indices(i) = 1;
        end
        
    % Is the length too big?
    elseif lengths(i) > cc_lengths(end)
        
        % Find the closest index.
        if (lengths(i) - cc_lengths(end)) / cc_lengths(1) < .5
            indices(i) = length(cc_lengths);
        else
            indices(i) = 1;
        end

    % Find the closest index.
    else
        
        % Try jumping to just before the requested length.
        % Note: most chain-code lengths advance by at most sqrt(2) at each
        % index. But I don't trust IEEE division so I use 1.5 instead.
        j = round(lengths(i) / 1.5) + 1;
        
        % Did we jump past the requested length?
        if j > length(cc_lengths) || lengths(i) < cc_lengths(j)
            j = 1;
        end
        
        % find the closest index.
        distJ = abs(lengths(i) - cc_lengths(j));
        while j < length(cc_lengths)
            
            % Is this index closer than the next one?
            % Note: overlapping points have equal distances. Therefore, if
            % the distances are equal, we advance.
            distNextJ = abs(lengths(i) - cc_lengths(j + 1));
            if distJ < distNextJ
                break;
            end
            
            % Advance.
            distJ = distNextJ;
            j = j + 1;
        end
        
        % Record the closest index.
        indices(i) = j;
    end
end

end

