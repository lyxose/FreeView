function normCurves = normalize_by_dim(subjCurves, normMode, normDim)
    %NORMALIZE_BY_dim Normalize subject curves along specified dimension.
    %
    %   normCurves = NORMALIZE_BY_dim(subjCurves, normMode, normDim)
    %
    %   输入参数:
    %       subjCurves - 多维数组，指定维度为时间或序列，其余维度会被展开为行。
    %       normMode   - 字符串，归一化方式，可选:
    %                       'zscore'  : 零均值单位方差归一化
    %                       'minmax'  : 最小-最大归一化到 [0,1]
    %                       'sum1'    : 总和归一化为 1
    %                       'demean'  : 去均值
    %                       'mean1'   : 均值归一化为 1
    %                       其他      : 不做归一化，直接返回原数据
    %       normDim    - (可选) 指定归一化的维度，默认最后一维
    %
    %   输出参数:
    %       normCurves - 与 subjCurves 形状相同的归一化结果。
    epsv = 1e-12;
    sz = size(subjCurves);
    if nargin < 3 || isempty(normDim)
        normDim = ndims(subjCurves);
    end
    nTime = sz(normDim);
    permOrder = [normDim, setdiff(1:numel(sz), normDim)];
    invOrder = zeros(1, numel(sz));
    invOrder(permOrder) = 1:numel(sz);
    subjMat = permute(subjCurves, permOrder);
    subjMat = reshape(subjMat, nTime, []);
    normMat = subjMat;
    for si = 1:size(normMat,2)
        r = normMat(:,si);
        switch lower(normMode)
            case 'zscore'
                mu = mean(r); sd = std(r);
                normMat(:,si) = (r - mu) ./ max(sd, epsv);
            case 'minmax'
                mn = min(r); mx = max(r);
                normMat(:,si) = (r - mn) ./ max(mx-mn, epsv);
            case 'sum1'
                s = sum(r);
                normMat(:,si) = r ./ max(s, epsv);
            case 'demean'
                mu = mean(r);
                normMat(:,si) = r - mu;
            case 'mean1'
                mu = mean(r);
                normMat(:,si) = r ./ max(mu, epsv);
            otherwise
                normMat(:,si) = r;
        end
    end
    normCurves = ipermute(reshape(normMat, [nTime, sz(setdiff(1:numel(sz), normDim))]), permOrder);
end
