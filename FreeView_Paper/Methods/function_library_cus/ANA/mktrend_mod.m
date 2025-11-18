
function [h, p, S, Z, sen, n_eff, acf1] = mktrend_mod(x, maxLag)
% 修正 Mann-Kendall 检验（Hamed & Rao, 1998），并计算 Sen's slope
% x: 输入序列
% maxLag: 计算自相关的最大滞后阶数（推荐10）

    x = x(:);
    x = x(~isnan(x));
    n = length(x);
    if nargin < 2, maxLag = round(4*(n/100)^(2/9)); end  % Newey–West自动带宽规则

    % 1. 计算 S
    S = 0;
    for k = 1:n-1
        S = S + sum(sign(x((k+1):n) - x(k)));
    end

    % 2. 计算原始方差
    unique_x = unique(x);
    g = length(unique_x);
    if n == g
        varS = n*(n-1)*(2*n+5)/18;
    else
        tp = histcounts(x, [unique_x; unique_x(end)+1]);
        varS = n*(n-1)*(2*n+5);
        for j = 1:length(tp)
            varS = varS - tp(j)*(tp(j)-1)*(2*tp(j)+5);
        end
        varS = varS / 18;
    end

    % 3. 计算自相关系数（ACF），排除lag=0
    acf = autocorr(x, min(maxLag, n-2)); % 需要Financial Toolbox，Econometrics Toolbox
    acf1 = acf(2:end); % lag=1:maxLag

    % 4. 计算有效样本量 n_eff
    sum_r = sum((n - (1:maxLag)) ./ (n * (n-1)) .* acf1);
    n_eff = n / (1 + 2*sum(sum_r)*n);

    % 5. 修正方差
    varS_mod = varS * n / n_eff;

    % 6. 计算Z
    if S > 0
        Z = (S - 1)/sqrt(varS_mod);
    elseif S < 0
        Z = (S + 1)/sqrt(varS_mod);
    else
        Z = 0;
    end

    % 7. p值
    p = 2 * (1 - normcdf(abs(Z), 0, 1));
    h = p < 0.05;

    % 8. Sen's slope
    slopes = [];
    for i = 1:n-1
        for j = i+1:n
            slopes(end+1) = (x(j) - x(i)) / (j - i);
        end
    end
    sen = median(slopes);

end
