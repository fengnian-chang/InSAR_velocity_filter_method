function B = loess_nanfit2(A, winsize, varargin)
% LOESS_NANFIT2: Weighted LOESS smoothing (local planar fit) with NaN support
%
% Syntax:
%   B = loess_nanfit2(A, winsize)
%   B = loess_nanfit2(A, winsize, 'nanout')
%   B = loess_nanfit2(A, winsize, 'nonanout')
%
% Inputs:
%   A        - Input 2D matrix (with possible NaNs)
%   winsize  - Odd integer window size (e.g., 11, 21, 51)
%   'nanout'   : Output retains NaNs in original A
%   'nonanout' : Output fills as many pixels as possible (default)
%
% Output:
%   B        - Smoothed output matrix

% Check winsize
if mod(winsize, 2) ~= 1
    error('winsize must be odd.');
end

% Parse mode
nanout = false;
if nargin > 2
    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'nanout')
            nanout = true;
        end
    end
end

% Initialization
[h, w] = size(A);
B = nan(h, w);
radius = floor(winsize / 2);

% Loop over pixels
for i = 1:h
    for j = 1:w
        % Define local window bounds
        r1 = max(1, i - radius); r2 = min(h, i + radius);
        c1 = max(1, j - radius); c2 = min(w, j + radius);

        % Extract window data
        W = A(r1:r2, c1:c2);
        [xgrid, ygrid] = meshgrid(c1:c2, r1:r2);
        x = xgrid(:); y = ygrid(:); z = W(:);

        % Keep valid points
        valid = ~isnan(z);
        x = x(valid); y = y(valid); z = z(valid);

        if numel(z) < 6
            continue;  % Not enough points
        end

        % Compute normalized distance to center (for weights)
        dx = x - j;
        dy = y - i;
        r = sqrt(dx.^2 + dy.^2);
        r_norm = r / max(radius, 1e-6);  % Normalize to [0,1]

        % Tricube weights
        wgt = (1 - r_norm.^3).^3;
        wgt(r_norm > 1) = 0;

        % Weighted least squares plane fit: z = a + bx + cy
        A_fit = [ones(size(x)), x, y];
        W = diag(wgt);
        if rank(A_fit) < 3
            continue;  % Ill-conditioned
        end

        coeff = (A_fit' * W * A_fit) \ (A_fit' * W * z);

        % Evaluate at center point (j, i)
        B(i,j) = [1, j, i] * coeff;
    end
end

% If 'nanout' mode, restore original NaNs
if nanout
    B(isnan(A)) = NaN;
end
end
