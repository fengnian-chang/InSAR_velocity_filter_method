function B = nanmedfilt2(A, winsize, mode)
% NANMEDFILT2 - NaN-aware 2D median filter with controllable NaN output
%
% Syntax:
%   B = nanmedfilt2(A, winsize)
%   B = nanmedfilt2(A, winsize, mode)
%
% Inputs:
%   A        - Input matrix with NaNs
%   winsize  - 2-element vector, e.g., [m n]
%   mode     - 'nanout' (default): preserve NaNs in output where A is NaN
%            - 'nonanout': fill all values (NaNs will be smoothed if possible)
%
% Output:
%   B - Filtered matrix

if nargin < 3
    mode = 'nanout';
end
assert(isequal(length(winsize),2), 'winsize must be a 2-element vector');

pad = floor(winsize/2);
A_pad = padarray(A, pad, NaN);
B = nan(size(A));

for i = 1:size(A,1)
    for j = 1:size(A,2)
        block = A_pad(i:i+2*pad(1), j:j+2*pad(2));
        valid = block(~isnan(block));
        if ~isempty(valid)
            B(i,j) = median(valid);
        end
    end
end

% Handle NaN mask based on mode
if strcmpi(mode, 'nanout')
    B(isnan(A)) = NaN;
end
end
