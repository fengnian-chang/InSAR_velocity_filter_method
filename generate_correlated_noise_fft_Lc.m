function noise2d = generate_correlated_noise_fft_Lc(X_km, Y_km, Lc_km, target_std)
% GENERATE_CORRELATED_NOISE_FFT_LC
% 使用 FFT 生成具有指定相关长度 Lc 的 2D 空间相关噪声。
% 这个版本允许显式控制平滑程度（Lc，单位 km）

    [ny, nx] = size(X_km);

    % 像元间距
    x_vec = X_km(1, :);
    y_vec = Y_km(:, 1);
    dx = mean(diff(x_vec));
    dy = mean(diff(y_vec));
    sp = mean([dx, dy]);   % km

    % 频率坐标（cycles per km）
    fx = (0:nx-1);
    fy = (0:ny-1);
    [FX, FY] = meshgrid(fx, fy);

    % shift to negative frequencies
    FX = ifftshift(FX - floor(nx/2));
    FY = ifftshift(FY - floor(ny/2));

    freq = sqrt(FX.^2 + FY.^2) ./ (nx * sp);   % cycles per km

    % ------------------------------
    % ★ 控制噪声平滑程度的关键：功率谱 S(f)
    %   采用 isotropic exponential model 的傅里叶对应形式：
    %
    %       S(f) ∝ 1 / (1 + (f * Lc)^2)
    %
    %   Lc 越大 → 高频衰减越大 → 更平滑
    % ------------------------------
    S = 1 ./ (1 + (freq * Lc_km).^2);

    % 对应振幅谱（sqrt(power spectrum)）
    A = sqrt(S);

    % 白噪声
    white_noise = randn(ny, nx);

    % 频域滤波
    Yf = fft2(white_noise);
    Yf_filt = Yf .* A;

    % 回到空间域
    y = real(ifft2(Yf_filt));

    % 标准化 & 缩放
    y = y - mean(y(:));
    y = y ./ std(y(:));
    noise2d = y * target_std;

end
