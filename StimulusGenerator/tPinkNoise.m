function truncatedNoise = tPinkNoise(N, seed, contrast, show)
    % 生成2D 1/f粉噪声纹理，灰度级，在2*sd处截断
    % Generate a 2D 1/f pink noise texture in gray level, truncated at 2*sd
    %
    % -----Parameters-----
    % N: 图像尺寸 (N x N)
    %    Image size (N x N)
    % seed: 随机种子
    %       Random seed
    % contrast: 对比度 (0-1), 默认1
    %           Contrast (0-1), default 1
    % show: 是否显示图像, 默认false
    %       Whether to display image, default false
    %
    % -----Returns-----
    % truncatedNoise: 粉噪声图像矩阵，范围[0, 1]
    %                 Pink noise image matrix, range [0, 1]
    %
    % Yuxin Lu, IPCAS
    % luyx@psych.ac.cn
    % 2025.1.7
    
    if nargin < 3
        contrast = 1;
    end
    
    if nargin < 4
        show = false;
    end
    
    currentSeed = rng;  % 保存当前随机种子状态 / Save current rng seed state
    rng(seed);

    % 创建频率网格
    % Create a grid of frequencies
    [u, v] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
    freq = sqrt(u.^2 + v.^2);  % 径向频率 / Radial frequency
    
    % 避免DC分量除零
    % Avoid division by zero for the DC component
    freq(freq == 0) = 1;
    
    % 生成随机相位
    % Generate random phases
    randomPhase = exp(1i * 2 * pi * rand(N, N));
    
    % 按1/f缩放幅度
    % Scale the amplitude by 1/f
    amplitude = 1 ./ freq;
    amplitude = amplitude / max(amplitude(:));  % 归一化幅度 / Normalize amplitude
    
    % 创建噪声的频域表示
    % Create the frequency-domain representation of the noise
    frequencyDomainNoise = amplitude .* randomPhase;
    
    % 转换回空间域
    % Transform back to the spatial domain
    spatialNoise = real(ifft2(ifftshift(frequencyDomainNoise)));
    
    % 计算标准差
    % Calculate standard deviation
    stdNoise = std(spatialNoise(:));
    
    % 在2*sd处截断
    % Truncate at 2*sd
    truncatedNoise = spatialNoise;
    truncatedNoise(abs(truncatedNoise) > 2*stdNoise) = ...
        2 * sign(truncatedNoise(abs(truncatedNoise) > 2*stdNoise)) * stdNoise;
    
    % 归一化到[0, 1]
    % Normalize to range [0, 1]
    truncatedNoise = (truncatedNoise - min(truncatedNoise(:))) / ...
                     (max(truncatedNoise(:)) - min(truncatedNoise(:)));
    
    % 应用对比度
    % Apply contrast
    truncatedNoise = contrast * (truncatedNoise - 0.5) + 0.5;
    
    % 显示1/f噪声图像
    % Display the 1/f noise image
    if show
        figure;
        imshow(truncatedNoise, []);
        title('1/f Pink Noise (1/f 粉噪声)');
        colorbar;
    end
    
    rng(currentSeed);  % 恢复之前的随机种子 / Reset rng as before

end
