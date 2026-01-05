function smoothed = gaussian_smooth_along_dim(mat, sigma, dim, mode)
    % 高斯平滑函数：对输入矩阵沿指定维度进行高斯平滑
    if nargin < 3, dim = 2; end
    if nargin < 4, mode = 'zero'; end
    radius = ceil(3 * sigma);
    tvec = -radius:radius;
    gk = exp(-0.5 * (tvec ./ sigma).^2);
    gk = gk / sum(gk);

    sz = size(mat);
    smoothed = zeros(sz, 'like', mat);

    idx = repmat({':'}, 1, ndims(mat));
    otherDims = setdiff(1:ndims(mat), dim);
    szOther = sz(otherDims);

    % Padding function
    function arr = pad_slice(slice)
        padsize = [0 0];
        padsize(dim) = radius;
        switch lower(mode)
            case 'circular'
                arr = padarray(slice, padsize, 'circular', 'both');
            case 'replicate'
                arr = padarray(slice, padsize, 'replicate', 'both');
            case 'symmetric'
                arr = padarray(slice, padsize, 'symmetric', 'both');
            otherwise % 'zero'
                arr = padarray(slice, padsize, 0, 'both');
        end
    end

    if isempty(otherDims)
        slice = mat;
        slice_pad = pad_slice(slice);
        smoothed = conv(slice_pad, gk, 'valid');
    else
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
            slice_pad = pad_slice(slice);
            smoothed(idx_sub{:}) = conv(slice_pad, gk, 'valid');
        end
    end
end