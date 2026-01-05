function B = nanmedfilt2(A, winsize, mode)
% NANMEDFILT2 - NaN-aware 2D median filter with rectangular or circular window
%
% Syntax:
%   B = nanmedfilt2(A, [m n])           % rectangular window, NaN-out (default)
%   B = nanmedfilt2(A, [m n], mode)     % rectangular window
%   B = nanmedfilt2(A, d)               % circular window with diameter d
%   B = nanmedfilt2(A, d, mode)         % circular window with diameter d
%
% Inputs:
%   A        - Input matrix with NaNs
%   winsize  - If scalar: circular window diameter (in pixels)
%              If 2-element vector [m n]: rectangular window size
%   mode     - 'nanout' (default): preserve NaNs in output where A is NaN
%            - 'nonanout': fill all values (NaNs will be smoothed if possible)
%
% Output:
%   B - Filtered matrix

if nargin < 3
    mode = 'nanout';
end

% Prepare output
B = nan(size(A));

% ---- Rectangular window (original behaviour): winsize is [m n] ----
if numel(winsize) == 2 && ~isscalar(winsize)

    winsize = double(winsize(:).');           % force row vector
    assert(numel(winsize) == 2, ...
        'winsize must be a scalar (for circular) or 2-element vector (for rectangular).');

    pad = floor(winsize/2);
    A_pad = padarray(A, pad, NaN);

    for i = 1:size(A,1)
        for j = 1:size(A,2)
            block = A_pad(i:i+2*pad(1), j:j+2*pad(2));
            valid = block(~isnan(block));
            if ~isempty(valid)
                B(i,j) = median(valid);
            end
        end
    end

% ---- Circular window: winsize is scalar diameter ----
elseif isscalar(winsize)

    % enforce integer, odd diameter
    d = round(double(winsize));
    assert(d >= 1, 'winsize (diameter) must be >= 1.');
    if mod(d,2) == 0
        d = d + 1;  % make it odd so there is a clear centre
    end

    r = (d - 1) / 2;               % radius in pixels
    pad = [r r];                   % symmetric padding
    A_pad = padarray(A, pad, NaN);

    % precompute circular mask
    [X, Y] = meshgrid(1:d, 1:d);
    cx = (d + 1) / 2;
    cy = (d + 1) / 2;
    circleMask = (X - cx).^2 + (Y - cy).^2 <= r^2 + eps;

    for i = 1:size(A,1)
        for j = 1:size(A,2)
            block = A_pad(i:i+2*pad(1), j:j+2*pad(2));
            valid = block(circleMask & ~isnan(block));
            if ~isempty(valid)
                B(i,j) = median(valid);
            end
        end
    end

else
    error('winsize must be a scalar (for circular) or 2-element vector (for rectangular).');
end

% Handle NaN mask based on mode
if strcmpi(mode, 'nanout')
    B(isnan(A)) = NaN;
end
end
