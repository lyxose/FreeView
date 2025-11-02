function smoothed = gaussian_smooth_along_dim(mat, sigma, dim)
    % 高斯平滑函数：对输入矩阵沿指定维度进行高斯平滑
    if nargin < 3, dim = 2; end
    radius = ceil(3 * sigma);
    tvec = -radius:radius;
    gk = exp(-0.5 * (tvec ./ sigma).^2);
    gk = gk / sum(gk);

    sz = size(mat);
    smoothed = zeros(sz, 'like', mat);

    % 按指定维度卷积
    idx = repmat({':'}, 1, ndims(mat));
    otherDims = setdiff(1:ndims(mat), dim);
    szOther = sz(otherDims);
    if isempty(otherDims)
        smoothed = conv(mat, gk, 'same');
    else
        % Generate all index combinations for otherDims (corrected)
        nOther = numel(otherDims);
        indexRanges = arrayfun(@(s) 1:s, szOther, 'UniformOutput', false);
        indexGrid = cell(1, nOther);
        [indexGrid{:}] = ndgrid(indexRanges{:});
        subsMat = cellfun(@(x) x(:), indexGrid, 'UniformOutput', false);
        nComb = numel(subsMat{1});
        for k = 1:nComb
            idx_sub = idx;
            for d = 1:nOther
                idx_sub{otherDims(d)} = subsMat{d}(k);
            end
            slice = squeeze(mat(idx_sub{:}));
            smoothed(idx_sub{:}) = conv(slice, gk, 'same');
        end
    end
end