function [c] = circConv(a, b, filter_length)
%circConv Convolve the circularly connected vector a with b.
%
%   c = circConv(a, b)
%
%   c = circConv(a, [], filter_length)
%
%   Inputs:
%       a - a circularly connected vector
%       b - the filter vector
%       filter_length - length of a box filter
%
%   Outputs:
%       c - the convolution of the circularly connected vector a with b
%

if isempty(b)
   b = 1/filter_length*ones(1,filter_length); 
end

% Are the inputs vectors.
if ~isvector(a) || ~isvector(b)
  error('circConv:AorBNotVector', 'A and B must be vectors');
end

% Is the convolution vector too long?
if length(a) < length(b)
  warning('circConv:AsmallerThanB', ...
      'A is smaller than B and, therefore, they cannot be convolved');
  c = a;
  return;
end

% Wrap the ends of A and convolve with B.
wrapSize = ceil(length(b) / 2);
% wrapA(1:wrapSize) = a((end - wrapSize + 1):end);
% wrapA((end + 1):(end + length(a))) = a;
% wrapA((end + 1):(end + wrapSize)) = a(1:wrapSize);

if size(a,2) > 1
    a = a';
end

wrapA = [a((end - wrapSize + 1):end); a; a(1:wrapSize)];
wrapB = conv(wrapA, b, 'same');

% Strip away the wrapped ends.
c = wrapB((wrapSize + 1):(end - wrapSize));
end

