function B = loess_nanfit2_plane(A, winsize, varargin)
% LOESS_NANFIT2: Apply LOESS (local planar fit) smoothing to 2D matrix with NaN support
%
%   B = loess_nanfit2_plane(A, winsize)
%   B = loess_nanfit2_plane(A, winsize, 'nanout')
%
% Inputs:
%   A        - Input 2D matrix, can contain NaN
%   winsize  - Size of moving window (must be odd integer)
%   'nanout' - (Optional) Keep NaN locations from input in output
%
% Output:
%   B        - Smoothed 2D output matrix (NaN-aware planar smoothing)
%
% Example:
%   B = loess_nanfit2_plane(V, 11, 'nanout');

% ---------------------
% Validate input
% ---------------------
if mod(winsize, 2) ~= 1
    error('winsize must be odd.');
end

nanout = false;
if nargin > 2
    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'nanout')
            nanout = true;
        end
    end
end

% ---------------------
% Initialize
% ---------------------
[h, w] = size(A);
B = nan(h, w);
radius = floor(winsize / 2);

% ---------------------
% Loop over pixels
% ---------------------
for i = 1:h
    for j = 1:w
        % Define window bounds
        r1 = max(1, i - radius); r2 = min(h, i + radius);
        c1 = max(1, j - radius); c2 = min(w, j + radius);

        % Extract local window and coordinates
        W = A(r1:r2, c1:c2);
        [xwin, ywin] = meshgrid(c1:c2, r1:r2);
        x = xwin(:); y = ywin(:); z = W(:);

        % Keep valid (non-NaN) points
        valid = ~isnan(z);
        x = x(valid); y = y(valid); z = z(valid);

        % Skip if insufficient data
        if numel(z) < 6
            continue;
        end

        % Fit plane: z = a + bx + cy
        A_fit = [ones(size(x)), x, y];
        if rank(A_fit) < 3
            continue; % Ill-conditioned
        end

        coeffs = A_fit \ z;

        % Evaluate fit at center point (j, i)
        B(i,j) = [1, j, i] * coeffs;
    end
end

% ---------------------
% Restore original NaN mask if needed
% ---------------------
if nanout
    B(isnan(A)) = NaN;
end
end
