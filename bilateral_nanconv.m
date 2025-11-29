function V_bilateral = bilateral_nanconv( ...
    V, spatial_sigma, range_sigma, ...
    spatial_winsize, center_winsize, ...
    center_mode, varargin)
% BILATERAL_NANCONV  Bilateral filter that handles NaN values in the input.
% 
%   V_bilateral = bilateral_nanconv(V, spatial_sigma, range_sigma, ...
%                                   spatial_winsize, center_winsize, ...
%                                   center_mode, ...)
%
%   center_mode 选项：
%     'median'   - range 的中心值采用局部中值（窗口大小 = center_winsize）
%                  调用示例：
%                     bilateral_nanconv(V, ss, rs, sw, cw, 'median', 'nanout')
%
%     'gaussian' - range 的中心值采用事先 Gaussian 平滑后的 V_center(i,j)
%                  Gaussian 核的窗口大小 = center_winsize，sigma = center_sigma
%                  调用示例：
%                     bilateral_nanconv(V, ss, rs, sw, cw, 'gaussian', center_sigma, 'nanout')
%
%   Optional flags（在 center_mode 参数之后）：
%     'nanout'   - 输出在输入为 NaN 的位置保持 NaN
%     'nonanout' - 能插值就插值，只有确实无法插值的地方才 NaN
%
%   INPUTS:
%     V               - Input 2D matrix with NaNs
%     spatial_sigma   - Std of spatial Gaussian kernel (in px)
%     range_sigma     - Std of range kernel (in data units, e.g. mm/yr)
%     spatial_winsize - Window size for spatial kernel (odd number)
%     center_winsize  - Window size for local center (median 或 Gaussian)
%     center_mode     - 'median' or 'gaussian'
%     varargin        - 对于 'gaussian'：第一个为 center_sigma，其余为 flags
%                       对于 'median'：直接就是 flags
%
%   OUTPUT:
%     V_bilateral     - Filtered result of same size as V

% ------------ 解析 center_mode & 额外参数 -------------
center_mode = lower(center_mode);

if strcmp(center_mode, 'gaussian')
    % 第一个额外参数必须是 center_sigma
    if isempty(varargin) || ~isnumeric(varargin{1})
        error('For center_mode="gaussian", first extra argument must be numeric center_sigma.');
    end
    center_sigma = varargin{1};
    flag_args    = varargin(2:end);
elseif strcmp(center_mode, 'median')
    % median 模式不需要 center_sigma
    center_sigma = [];
    flag_args    = varargin;
else
    error('center_mode must be "median" or "gaussian".');
end

% ------------ 解析 flags: nanout / nonanout -------------
nanout = false;
for k = 1:numel(flag_args)
    switch lower(flag_args{k})
        case 'nanout'
            nanout = true;
        case 'nonanout'
            nanout = false;
    end
end

% ------------ 初始化 -------------
[m, n] = size(V);
V_bilateral = nan(m, n);

valid_mask_global = ~isnan(V);

% ------------ 空间核（Gaussian for spatial weights）-------------
half_win_spatial = floor(spatial_winsize/2);
[Xg, Yg] = meshgrid(-half_win_spatial:half_win_spatial, -half_win_spatial:half_win_spatial);
spatial_kernel = exp(-(Xg.^2 + Yg.^2) / (2 * spatial_sigma^2));

% ------------ 若为 gaussian 模式，预先计算 V_center -------------
if strcmp(center_mode, 'gaussian')
    half_win_center = floor(center_winsize/2);
    [Xc, Yc] = meshgrid(-half_win_center:half_win_center, -half_win_center:half_win_center);
    center_kernel = exp(-(Xc.^2 + Yc.^2) / (2 * center_sigma^2));
    center_kernel = center_kernel / sum(center_kernel(:));
    
    % nan-aware 卷积得到 V_center
    V_center = nanconv(V, center_kernel, 'edge', 'nanout');
end

% median 模式下的中心窗口半径
if strcmp(center_mode, 'median')
    half_win_center = floor(center_winsize/2);
end

% ------------ 主循环 -------------
for i = 1:m
    for j = 1:n
        
        % === 1. 空间窗口 ===
        r1 = max(i - half_win_spatial, 1);
        r2 = min(i + half_win_spatial, m);
        c1 = max(j - half_win_spatial, 1);
        c2 = min(j + half_win_spatial, n);
        
        local_patch = V(r1:r2, c1:c2);
        local_valid = valid_mask_global(r1:r2, c1:c2);
        
        % 对应截取 spatial kernel
        sk = spatial_kernel( ...
            (r1 - i + half_win_spatial + 1):(r2 - i + half_win_spatial + 1), ...
            (c1 - j + half_win_spatial + 1):(c2 - j + half_win_spatial + 1) );
        
        % === 2. 定义 range 的中心值 center_val ===
        if strcmp(center_mode, 'gaussian')
            % 直接取预先算好的 Gaussian 平滑场
            center_val = V_center(i,j);
            
            if isnan(center_val)
                if nanout
                    V_bilateral(i,j) = NaN;
                    continue;
                else
                    V_bilateral(i,j) = NaN;
                    continue;
                end
            end
            
        elseif strcmp(center_mode, 'median')
            % 使用局部窗口的中值
            rr1 = max(i - half_win_center, 1);
            rr2 = min(i + half_win_center, m);
            cc1 = max(j - half_win_center, 1);
            cc2 = min(j + half_win_center, n);
            
            median_patch = V(rr1:rr2, cc1:cc2);
            median_valid = ~isnan(median_patch);
            
            if any(median_valid(:))
                center_val = median(median_patch(median_valid));
            else
                if nanout
                    V_bilateral(i,j) = NaN;
                    continue;
                else
                    V_bilateral(i,j) = NaN;
                    continue;
                end
            end
        end
        
        % === 3. range 权重，相对 center_val ===
        range_weight = zeros(size(local_patch));
        range_weight(local_valid) = exp( -((local_patch(local_valid) - center_val).^2) / (2 * range_sigma^2) );
        
        % 总权重 = 空间核 * range 核，并对 NaN 置零
        total_weight = sk .* range_weight;
        total_weight(~local_valid) = 0;
        
        denom = sum(total_weight(:));
        if denom > 0
            V_bilateral(i,j) = sum(local_patch(local_valid) .* total_weight(local_valid)) / denom;
        else
            V_bilateral(i,j) = NaN;
        end
        
        % === 4. nanout 逻辑：输入 NaN 的位置强制保持 NaN ===
        if nanout && isnan(V(i,j))
            V_bilateral(i,j) = NaN;
        end
    end
end

end
