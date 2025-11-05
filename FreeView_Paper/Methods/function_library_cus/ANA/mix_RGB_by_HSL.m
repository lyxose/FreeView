function mixed_rgb = mix_RGB_by_HSL(rgb1, rgb2)
    % 输入颜色为RGB格式，范围为 [0, 1]
    % 输出混合后的RGB颜色，范围为 [0, 1]

    % 转换为 [0, 1] 范围
    rgb1 = double(rgb1);
    rgb2 = double(rgb2);

    % 转换为 HSL
    hsl1 = rgb2hsl(rgb1);
    hsl2 = rgb2hsl(rgb2);

    % 混合色相（Hue）使用向量平均法
    h1 = hsl1(1) * 2 * pi;
    h2 = hsl2(1) * 2 * pi;
    x = cos(h1) + cos(h2);
    y = sin(h1) + sin(h2);
    avg_hue = atan2(y, x) / (2 * pi);
    if avg_hue < 0
        avg_hue = avg_hue + 1;
    end

    % 混合饱和度和亮度
    avg_sat = (hsl1(2) + hsl2(2)) / 2;
    avg_lum = (hsl1(3) + hsl2(3)) / 2;

    % 合并为混合HSL
    mixed_hsl = [avg_hue, avg_sat, avg_lum];

    % 转换回RGB
    mixed_rgb = hsl2rgb(mixed_hsl);
end

function hsl = rgb2hsl(rgb)
    r = rgb(1); g = rgb(2); b = rgb(3);
    maxc = max([r, g, b]);
    minc = min([r, g, b]);
    l = (maxc + minc) / 2;

    if maxc == minc
        h = 0; s = 0;
    else
        d = maxc - minc;
        s = d / (1 - abs(2 * l - 1));
        if maxc == r
            h = mod((g - b) / d, 6);
        elseif maxc == g
            h = ((b - r) / d) + 2;
        else
            h = ((r - g) / d) + 4;
        end
        h = h / 6;
    end
    hsl = [h, s, l];
end

function rgb = hsl2rgb(hsl)
    h = hsl(1); s = hsl(2); l = hsl(3);
    c = (1 - abs(2 * l - 1)) * s;
    x = c * (1 - abs(mod(h * 6, 2) - 1));
    m = l - c / 2;

    if h < 1/6
        r1 = c; g1 = x; b1 = 0;
    elseif h < 2/6
        r1 = x; g1 = c; b1 = 0;
    elseif h < 3/6
        r1 = 0; g1 = c; b1 = x;
    elseif h < 4/6
        r1 = 0; g1 = x; b1 = c;
    elseif h < 5/6
        r1 = x; g1 = 0; b1 = c;
    else
        r1 = c; g1 = 0; b1 = x;
    end

    rgb = [r1 + m, g1 + m, b1 + m];
end
