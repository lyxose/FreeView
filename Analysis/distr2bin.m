function bin_values = distr2bin(distr, img_width, img_height, edges_x, edges_y)
[x,y]=meshgrid(1:img_width, 1:img_height);
bin_values = zeros(length(edges_x)-1, length(edges_y)-1);
for i = 1:(length(edges_x)-1)
    for j = 1:(length(edges_y)-1)
        % 找到当前 bin 内的元素
        bin_mask = (x >= edges_x(i)) & (x < edges_x(i + 1)) & ...
                   (y >= edges_y(j)) & (y < edges_y(j + 1));
        % 累加当前 bin 内的元素值
        bin_values(i, j) = sum(distr(bin_mask));
    end
end
bin_values = bin_values';
end