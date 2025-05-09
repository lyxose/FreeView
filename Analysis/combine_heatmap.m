% this demo code is part of Titta, a toolbox providing convenient access to
% eye tracking functionality using Tobii eye trackers
%
% Titta can be found at https://github.com/dcnieho/Titta. Check there for
% the latest version.
% When using Titta, please cite the following paper:
%
% Niehorster, D.C., Andersson, R. & Nystrom, M., (2020). Titta: A toolbox
% for creating Psychtoolbox and Psychopy experiments with Tobii eye
% trackers. Behavior Research Methods.
% doi: https://doi.org/10.3758/s13428-020-01358-8

clear variables; clear global; clear mex; close all; fclose('all'); clc

dbstop if error % for debugging: trigger a debug point when an error occurs

% setup directories
myDir = fileparts(mfilename('fullpath'));
cd(myDir);
cd ..;                                  dirs.home               = cd;
    cd Analysis;                            dirs.ana            = cd;
            cd Processed_data;                  dirs.data       = cd;
                    cd samples_ophak;           dirs.samples    = cd;
            cd ..;  cd fixDet;                  dirs.fix        = cd;
            cd ..;  cd msgs_ophak;              dirs.msgsO      = cd;
            cd ..;  
    cd ..;  cd function_library;        dirs.funclib            = cd;  
    cd ..;  cd function_library_cus;    dirs.funclib_cus        = cd;
    cd ..;
cd ..;
    cd Data; 
            cd Formal;                  dirs.mat                = cd;
    cd ..;     
cd(dirs.ana);
addpath(genpath(dirs.funclib));                 % add dirs to path
addpath(genpath(dirs.funclib_cus));

binSize = 25; % in pixel
dotsize = 10;

skip_corr=false;% ignore correct trials
drawTraj = true;
drawEach = true;
    tagNum = true;
drawCorrDist = true;
fitModel = false;
% only plot the selected fixations
Start_nFix  = 1;
End_nFix    = 2; 
atLeast_nFix = End_nFix+5;
% for Start_nFix = 2:3
% for End_nFix = Start_nFix:3

R_max = 9.5;
R_min = 0;

%%% get result table with eye data
[files,~] = FileFromFolder(dirs.fix,[],'mat');
filtstr = '^Dat_Sub(\d+)_Ses(\d+).mat$';
matched = regexpi({files.name}.',filtstr,'tokens');
resfile_idx = ~cellfun(@isempty,matched);
resfiles   = files(resfile_idx);
sub_ses_res = cellfun(@(x) [str2double(x{1}{1}), str2double(x{1}{2})], matched(resfile_idx), 'UniformOutput', false);
sub_ses_res = cell2mat(sub_ses_res(~cellfun(@isempty, matched(resfile_idx))));

% select_sess = 1:length(resfiles);
select_sess = 2;
exclude_sess = [];
% exclude_sess = [1,2,3,10,18,25];
% select_sess = 10:14;
% select_sess = [1:10 14:length(resfiles)]; % single subj
% select_sess = find((sub_ses_res(:,2)==1)==1)';
% select_sess = [2,3,4,5,7,8,10,11,16,17];

alpha = 0.05; % traj

figure(1);  % 在循环外创建 figure，固定窗口
% figure(2);  % 在循环外创建 figure，固定窗口

% heatpath = sprintf('./results/comb_heatmap-%.0f-%.0fto%.0f/',length(select_sess),Start_nFix, End_nFix);
heatpath = sprintf('./results/comb_heatmap/');
trajpath = sprintf('./results/comb_trajec-%.0f-%.0fto%.0f/',length(select_sess),Start_nFix, End_nFix);
corrpath = sprintf('./results/comb_corr-%.0f-%.0fto%.0f/',length(select_sess),Start_nFix, End_nFix);
distpath = sprintf('./results/comb_dist-%.0f-%.0fto%.0f/',length(select_sess),Start_nFix, End_nFix);

if ~exist(heatpath,'dir')
    mkdir(heatpath);
end
% if ~exist(trajpath,'dir')
%     mkdir(trajpath);
% end
% if ~exist(corrpath,'dir')
%     mkdir(corrpath);
% end
% if ~exist(distpath,'dir')
%     mkdir(distpath);
% end
fixPos = [];
stTs = [];   % color will be depended on the start time 
triSpl = []; % to split each trial by the end index
tgLoc = [];
% 
% for ecc = [2,4,6]
% for ori = 0:45:315
%% 提取和合并数据
for p=select_sess
    if ismember(p,exclude_sess)
        continue
    end
    sess = load(fullfile(dirs.mat,sprintf("Dat_Sub%.0f_Ses%.0f.mat",sub_ses_res(p,1),sub_ses_res(p,2))),'expt','geometry');
    resT = load(fullfile(dirs.fix,[resfiles(p).fname '.mat'])).expT;
    img_width = sess.expt.winRect(3);  % 二维空间的宽度
    img_height = sess.expt.winRect(4); % 二维空间的高度
    ut = UT(sess.geometry.displayArea.width/10, img_width, mean(resT.headDist), false);
% for ori
%         rows = all(bsxfun(@eq, [resT.ECC, resT.Orient], [ecc,ori]), 2); 
%     for k = 1:max(ic)
%         fixPos = [];
%         rows = ~isnan(resT.judge);
%         iresT = resT(rows,:); % subset of resT, all in the EO ECC-Ori condition 
        iresT = resT;
        for i=1:height(iresT)
            if skip_corr && ~isnan(iresT.key2RT(i))
                continue
            end
            tFixPos = transpose([iresT.dat(i).fix.xpos; iresT.dat(i).fix.ypos]); % this trial eye trajectory
            tTime = iresT.dat(i).fix.startT;  % start time stamp
            tFixPos = tFixPos(tTime>0,:);
            tTime = tTime(tTime>0);
            % 剔除一开始就在中心1°以内的trial，并跳过长度等于1的trial
%             while length(tFixPos(:,1))>1 && norm(tFixPos(1,:)-[img_width, img_height]/2)<ut.deg2pix(1)
%                 if length(tFixPos(:,1))==1
%                     break
%                 end
%                 tFixPos = tFixPos(2:end,:); % remove the start fixation
%                 tTime = tTime(2:end,:);
%             end
            
            % cut to early fixations
            % 按剔除后总fixation数量筛选trial，以免选取到最后几个已经位于目标附近的trial
            if length(tTime(:,1))>=atLeast_nFix  %
                % exclude some trial that 被试分心 or 眨眼 导致fixation信号开始晚于1s
                if tTime(1)>1000
    %                 sprintf("Dat_Sub%.0f_Ses%.0f.mat",sub_ses_res(p,1),sub_ses_res(p,2))
    %                 sprintf('i=%.0f',i)
    %                 disp(iresT.dat(i).fix.startT);
                    continue
                end
                if End_nFix == 0
                    tFixPos = tFixPos(Start_nFix:end,:);
                    tTime = tTime(Start_nFix:end,:);
                else
%                     tFixPos = tFixPos(Start_nFix+7:End_nFix+7,:)-tFixPos(Start_nFix:End_nFix,:)+ones(End_nFix-Start_nFix+1,2).*[img_width/2,img_height/2];                
                    tFixPos = tFixPos(Start_nFix:End_nFix,:);
                    tTime = tTime(Start_nFix:End_nFix,:);
                end
            else
                continue % drop short trials
            end

%             tFixPos = arrayfun(@(i) ut.Pol2Rect([1,-1].*(ut.Rect2Pol(tFixPos(i, :)-[img_width, img_height]/2))-[0,ori]).*[1,-1]+[img_width, img_height]/2, ...
%                         1:size(tFixPos, 1), 'UniformOutput', false);
%             tFixPos = cell2mat(tFixPos');
            % concate all the selected results
            tsp(1)  = length(stTs)+1;
            stTs = [stTs;tTime];
            fixPos = [fixPos; tFixPos];
            tsp(2) = length(stTs);  
            triSpl = [triSpl;tsp];
            ttgLoc = ut.deg2pix([iresT.Xtarg(i),iresT.Ytarg(i)]).*[1,-1]+[img_width,img_height]./2;
            tgLoc = [tgLoc;ttgLoc];
        end
end
% end
% end

% figTitle = 
%% 绘制无序热图
        figure(1);
        clf;  
        edges_x = linspace(1, img_width,  round(img_width/binSize));  % 网格的X边界
        edges_y = linspace(1, img_height, round(img_height/binSize)); % 网格的Y边界
        counts = hist3(fixPos, 'Edges', {edges_x(1:end-1), edges_y(1:end-1)});

%         % statistics for the saccade distance
%         disPos = fixPos(2:end,:)-fixPos(1:end-1,:);
%         disPos(triSpl(:,2),:)=NaN;
%         disPos=disPos(~isnan(disPos(:,1)),:)+[img_width/2,img_height/2];
%         counts = hist3(disPos, 'Edges', {edges_x(1:end-1), edges_y(1:end-1)});
        density_map = counts' / sum(counts(:));
%         figure;
        imagesc(edges_x, edges_y, density_map);  % 使用imagesc绘制热图
        pr = prctile(density_map(:), 99.99);
        clim([min(density_map(:)), pr]);
%         set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        axis image;
        colorbar; % 显示颜色条
        % 标记中央位置：用 rectangle 函数绘制一个圆形
        hold on;  % 保持当前图形
        rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
                  'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 
%         rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
%                   'EdgeColor', 'r', 'LineWidth', 0.7, 'Curvature', [1, 1]);  %        
        % 标记目标位置
        for ecc = [2,4,6]
        for ori = 0:45:315
            tgWidth = ut.deg2pix(sess.expt.GaborWidth);
            target_loc = ut.Pol2Rect([ecc,ori]);
            target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
    
            rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
                      'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
        end
        end
        % 标记target范围
        tgECC_max = ut.deg2pix(R_max);
%         tgECC_min = ut.deg2pix(R_min);
%         rectangle('Position', [img_width/2-tgECC_min, img_height/2-tgECC_min, 2*tgECC_min, 2*tgECC_min], ...
%                   'EdgeColor', 'w', 'LineWidth', 0.7, 'Curvature', [1, 1]);
        rectangle('Position', [img_width/2-tgECC_max, img_height/2-tgECC_max, 2*tgECC_max, 2*tgECC_max], ...
                  'EdgeColor', 'w', 'LineWidth', 0.7, 'Curvature', [1, 1]);

        % 标记背景范围
        bgWidth = ut.deg2pix(sess.expt.bgWidth);
        rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
                  'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
        title(sprintf('Fixation heatmap %.0f to %.0f  %.0f subj'             ,Start_nFix,End_nFix,length(select_sess)-length(exclude_sess)));
        exportgraphics(gcf, sprintf([heatpath,'%.0fto%.0f_%.0fsess_heat.png'],Start_nFix,End_nFix,length(select_sess)-length(exclude_sess)), 'Resolution', 600);
        savefig(sprintf([heatpath,'%.0fto%.0f_%.0fsess_heat.fig']            ,Start_nFix,End_nFix,length(select_sess)-length(exclude_sess)))
        hold off; 
%%
        if drawTraj
            figure(2); clf; hold on; axis equal; 
            colormap(cool);  % 选择traj的颜色映射，例如 'parula', 'jet', 'hot'
            clim([0,1]);
            axis([0 img_width 0 img_height]);  % 横纵坐标范围为 img_width 和 img_height
            set(gca, 'YDir', 'reverse');  % y轴零点在左上角
            
            % 归一化时间值到 [0,1]，用于颜色映射
    %             normTime = time / maxTime;
    %         normTime = tTime / max(tTime);
            rankTime = tiedrank(stTs);
            normTime = rankTime/max(rankTime);
            
            % 绘制渐变线段/点
    %         for i = 1:size(triSpl,1)
    %             tsp = triSpl(i,:);
    %             xpos = fixPos(tsp(1):tsp(2),1);
    %             ypos = fixPos(tsp(1):tsp(2),2);
    %             tTime = normTime(tsp(1):tsp(2));
    %             for j = 1:length(xpos)%-1
    %                 ct = tTime(j);  % 当前时间的归一化值
    %                 color = cool(256);  % 获取 colormap 颜色
    %                 idx = max(1, round(ct * 255));  % 找到颜色索引
    %                 scatter(xpos(j), ypos(j), dotsize, color(idx, :), 'filled', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
    % %                 plot(xpos(j:j+1), ypos(j:j+1), '-', 'Color', [color(idx, :), alpha], 'LineWidth', 1.5);
    %             end
    %         end
            % 预计算所有点的颜色
            color = cool(256);  % 获取 colormap 颜色
            all_xpos = [];
            all_ypos = [];
            all_colors = [];
            
            for i = 1:size(triSpl, 1)
                tsp = triSpl(i, :);
                xpos = fixPos(tsp(1):tsp(2), 1);
                ypos = fixPos(tsp(1):tsp(2), 2);
                tTime = normTime(tsp(1):tsp(2));
                
                % 预计算每个点的颜色
                idx = max(1, round(tTime * 255));  % 找到颜色索引
                colors = color(idx, :);
                % 逐个检查轨迹
                if drawEach
                    if tagNum
                        for ii = 1:length(xpos)
                            text(xpos(ii), ypos(ii), num2str(ii), ...
                                 'FontSize', 5, 'Color', colors(ii,:));
                        end
                    else
                        scatter(xpos, ypos, dotsize, colors, 'filled');
                    end
                    rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
                              'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
                    rectangle('Position', [tgLoc(i,:)-5, 10, 10], ...
                              'EdgeColor', 'r', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
                    
                    input(sprintf('Next:%d',i+1))
                    cla;
                end
                
                % 收集所有点的坐标和颜色
                all_xpos = [all_xpos; xpos];
                all_ypos = [all_ypos; ypos];
                all_colors = [all_colors; colors];
            end
            
            % 一次性绘制所有点
            scatter(all_xpos, all_ypos, dotsize, all_colors, 'filled', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
    
            % 设置图像的坐标范围和 y 轴的零点在左上角
            % 标记中心位置
            rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
                      'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
            % 标记目标位置
            for ecc = [2,4,6]
            for ori = 0:45:315
                tgWidth = ut.deg2pix(sess.expt.GaborWidth);
                target_loc = ut.Pol2Rect([ecc,ori]);
                target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
        
                rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
                          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);
            end
            end
            % 标记背景范围
            bgWidth = ut.deg2pix(sess.expt.bgWidth);
            rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
                      'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    
            c = colorbar;  % 显示颜色条
            qtls = prctile(stTs,round(c.Ticks*100));
            c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), qtls, 'UniformOutput', false);
            c.Label.String = 'Start Time (ms)';
            title(sprintf('Fixation Trajectory %.0f to %.0f  %.0f subj',Start_nFix,End_nFix,length(select_sess)));
    %         xlabel('X Position');
    %         ylabel('Y Position');
            exportgraphics(gcf, sprintf([trajpath,'%.0fto%.0f_%.0fsess_tra.png'],Start_nFix,End_nFix,length(select_sess)), 'Resolution', 600);
            savefig(sprintf([trajpath,'%.0fto%.0f_%.0fsess_tra.fig'],Start_nFix,End_nFix,length(select_sess)))
            hold off; 
        end
        
        if drawCorrDist
            % bin的中心坐标
            [ex,ey]=meshgrid((edges_x(1:end-1)+edges_x(2:end))/2,...
                             (edges_y(1:end-1)+edges_y(2:end))/2);
            ex = ex-img_width/2;
            ey = ey-img_height/2;
            ec = sqrt(ex.^2+ey.^2);
            cir_mask = ec<bgWidth/2 & ec>ut.deg2pix(1);
            obs_data=density_map(cir_mask);
            obs_data=obs_data./sum(obs_data(:));
            exp_data_uni = ones(length(obs_data),1)./length(obs_data);
            sigma_space = 1.5.^(-10:25);
            corr_4gau = zeros(length(sigma_space),1);
            dist_4gau = zeros(length(sigma_space),1);
            corr_2dgau = zeros(length(sigma_space),1);
            dist_2dgau = zeros(length(sigma_space),1);
            
            if fitModel          
            h = waitbar(0, 'Processing...'); % 初始进度为 0            
            for k = 1:length(sigma_space)
            
            sigma = sigma_space(k);
            hold on
%             exp_distr_4gau = gaumap(img_width, img_height, sigma);
            sigma_r = 200;
            exp_distr_4gau = OctagonalPolmap(img_width, img_height, sigma_r, sigma);
            bin_values_4gau = distr2bin(exp_distr_4gau, img_width, img_height, edges_x, edges_y);
            bin_values_4gau(~cir_mask) = 0;
            bin_values_4gau = bin_values_4gau./sum(bin_values_4gau(:));
            exp_data_4gau = bin_values_4gau(cir_mask);
            corr_4gau(k) = corr(exp_data_4gau,obs_data);
            dist_4gau(k) = norm(exp_data_4gau-obs_data);
            
            exp_distr_2dgau=  fspecial('gaussian', [img_height, img_width], sigma);
            bin_values_2dgau = distr2bin(exp_distr_2dgau, img_width, img_height, edges_x, edges_y);
            bin_values_2dgau(~cir_mask) = 0;
            bin_values_2dgau = bin_values_2dgau./sum(bin_values_2dgau(:));
            exp_data_2dgau = bin_values_2dgau(cir_mask);
            corr_2dgau(k) = corr(exp_data_2dgau,obs_data);
            dist_2dgau(k) = norm(exp_data_2dgau-obs_data);
            
            % imagesc(edges_x, edges_y, bin_values_4gau);
            % rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
            %           'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
            % rectangle('Position', [img_width/2-5, img_height/2-5, 10, 10], ...
            %           'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
            % axis image;
            % colorbar; % 显示颜色条
            % hold off
            waitbar(k / length(sigma_space), h, sprintf('Processing... %d%%', round(k / length(sigma_space) * 100)));
            end
            
            figure(4); clf; hold on; 
            plot(sigma_space,corr_4gau , 'LineWidth', 1.5)
            plot(sigma_space,corr_2dgau, 'LineWidth', 1.5)
            plot(sigma_space,ones(k,1).*corr(exp_data_uni,obs_data), 'LineWidth', 1.5)
            legend('4gau','2dgau','uniform')
            % 设置横坐标为对数刻度
            set(gca, 'XScale', 'log');
            xlabel('sigma')
            ylabel('Correlation')
            title(sprintf('Heatmap correlation %.0f to %.0f  %.0f subj',Start_nFix,End_nFix,length(select_sess)));
            exportgraphics(gcf, sprintf([corrpath,'%.0fto%.0f_%.0fsess_corr.png'],Start_nFix,End_nFix,length(select_sess)), 'Resolution', 600);
            hold off
            
            
            figure(5); clf; hold on; 
            plot(sigma_space,dist_4gau , 'LineWidth', 1.5)
            plot(sigma_space,dist_2dgau, 'LineWidth', 1.5)
            plot(sigma_space,ones(k,1).*norm(exp_data_uni-obs_data), 'LineWidth', 1.5)
            legend('4gau','2dgau','uniform')
            % 设置横坐标为对数刻度
            set(gca, 'XScale', 'log');
            xlabel('sigma')
            ylabel('Euclidean distance')
            title(sprintf('Heatmap distance %.0f to %.0f  %.0f subj',Start_nFix,End_nFix,length(select_sess)));
            exportgraphics(gcf, sprintf([distpath,'%.0fto%.0f_%.0fsess_dist.png'],Start_nFix,End_nFix,length(select_sess)), 'Resolution', 600);
            hold off
            end

        end
 
% end
% end
% end
% end
if fitModel
    sigma=200;
    exp_distr_4gau = OctagonalPolmap(img_width, img_height, 100, 0.265, 100);
    bin_values_4gau = distr2bin(exp_distr_4gau, img_width, img_height, edges_x, edges_y);
    bin_values_4gau(~cir_mask) = 0;
    bin_values_4gau = bin_values_4gau./sum(bin_values_4gau(:));
    exp_data_4gau = bin_values_4gau(cir_mask);
    corr_4gau_best = corr(exp_data_4gau,obs_data)
    dist_4gau_best = norm(exp_data_4gau-obs_data)
    figure(11);
    imagesc(edges_x, edges_y, bin_values_4gau);
    rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
              'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    rectangle('Position', [img_width/2-5, img_height/2-5, 10, 10], ...
              'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    % 标记目标位置
    for ecc = [2,4,6]
    for ori = 0:45:315
        tgWidth = ut.deg2pix(sess.expt.GaborWidth);
        target_loc = ut.Pol2Rect([ecc,ori]);
        target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
        rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
                  'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    end
    end
    axis image;
    title('Σ4 Gau')
    colorbar; % 显示颜色条
    
    corr_4gau(k) = corr(exp_data_4gau,obs_data);
    
    exp_distr_2dgau=  fspecial('gaussian', [img_height, img_width], sigma);
    bin_values_2dgau = distr2bin(exp_distr_2dgau, img_width, img_height, edges_x, edges_y);
    bin_values_2dgau(~cir_mask) = 0;
    bin_values_2dgau = bin_values_2dgau./sum(bin_values_2dgau(:));
    exp_data_2dgau = bin_values_2dgau(cir_mask);
    corr_2dgau(k) = corr(exp_data_2dgau,obs_data);
    dist_2dgau(k) = norm(exp_data_2dgau-obs_data);
    figure(7);
    imagesc(edges_x, edges_y, bin_values_2dgau);
    rectangle('Position', [img_width/2-bgWidth/2, img_height/2-bgWidth/2, bgWidth, bgWidth], ...
              'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    rectangle('Position', [img_width/2-5, img_height/2-5, 10, 10], ...
              'EdgeColor', '#E3170D', 'LineWidth', 0.7, 'Curvature', [1, 1]);
    axis image;
    title('2d Gau')
    colorbar; % 显示颜色条
end

%     c == [ECC, Ori]
%     
%     % 假设dat.fix.xpos和dat.fix.ypos是1x51的double类型数据
%     % 设定图像的尺寸和分辨率
%     img_width = expt.winRect(3);  % 二维空间的宽度
%     img_height = expt.winRect(4); % 二维空间的高度
%     
%     % 假设dat.fix.xpos和dat.fix.ypos是1x51的double类型数据
%     xpos = dat.fix.xpos;
%     ypos = dat.fix.ypos;
%     % 1. 对二维平面中的点进行核密度估计（KDE），以估算每个位置的点出现概率
%     % 创建一个网格来表示二维空间
%     edges_x = linspace(1, img_width, 50);  % 网格的X边界
%     edges_y = linspace(1, img_height, 50); % 网格的Y边界
%     
%     % 2. 使用二维直方图，计算数据在网格上的分布
%     counts = hist3([xpos', ypos'], 'Edges', {edges_x, edges_y});
%     
%     % 3. 将结果转换为二维密度
%     density_map = counts' / sum(counts(:));  % 将计数转换为密度
%     
%     % 4. 绘制热图
%     figure;
%     imagesc(edges_x, edges_y, density_map);  % 使用imagesc绘制热图
%     axis xy; % 使y轴方向与图像坐标系一致
%     colorbar; % 显示颜色条
%     title('Density Heatmap');
%     xlabel('X Position');
%     ylabel('Y Position');
% 
% end
% 
% 
% 
% 
% 
% 
% 
% sessionFileName = sprintf('%s.mat',files(1).subj);
% expt = load(fullfile(dirs.mat,sessionFileName),'expt').expt;
% geom = load(fullfile(dirs.mat,sessionFileName),'geometry').geometry.displayArea;
% if nfiles~=height(res_Table)
%     error('some trials are losted')
% end
% 
% for p=1:nfiles
%     % load fix data
%     dat  = load(fullfile(dirs.fix,[files(p).fname '.mat'])); dat = dat.dat;
%     if isempty(dat.time)
%         warning('no data for %s, empty file',files(p).fname);
%         continue;
%     end
%     
%     % get msgs
%     msgs    = loadMsgs(fullfile(dirs.msgsO,[files(p).fname '.txt']));
%     [times,what,msgs] = parseMsgs_Form(msgs);
% 
% %     if ~strcmp(lastRead,sessionFileName)
% %         lastRead = sessionFileName;
% %         fInfo = [sess.expt.stim.fInfo];
% %     end
% %     qWhich= strcmp({fInfo.name},what{1});
%     
%     % load img, if only one
% %     if ~~exist(fullfile(dirs.AOImasks,what{1}),'file')
% %         img.data = imread(fullfile(dirs.AOImasks,what{1}));
% %     elseif ~~exist(fullfile(sess.expt.stim(qWhich).fInfo.folder,what{1}),'file')
% %         img.data = imread(fullfile(sess.expt.stim(qWhich).fInfo.folder,what{1}));
% %     else
% %     grayImg = uint8(sess.expt.stim*255);
%     ut = UT(geom.width/10, expt.winRect(3), res_Table.headDist(p));
%     bgCenter = expt.winRect(3:4)/2;
%     tgCenter_ = [res_Table.ECC(p),res_Table.Orient(p)];
%     tgCenter = ut.Pol2Rect(tgCenter_);
%     
%     grayImg = genStim(expt.winRect, ut, res_Table.bgContrast(p), ...
%         1, tgCenter, expt.GaborSF, expt.GaborWidth, ...
%         expt.GaborOrient, expt.bgWidth, res_Table.seed(p));
%     img.data     = uint8 (cat(3,grayImg,grayImg,grayImg)*255);
% 
% 
% 
% %     if ~isempty(img)
% %         % get position on screen
% %         stimRect = [0,0,size(img.data,2),size(img.data,1)];
% %         img.x    = linspace(stimRect(1),stimRect(3),size(img.data,2));
% %         img.y    = linspace(stimRect(2),stimRect(4),size(img.data,1));
% %     end
% %     
% %     % plot
% %     if ~ishghandle(fhndl)
% %         fhndl = figure('Units','normalized','Position',[0 0 1 1]);  % make fullscreen figure
% %     else
% %         figure(fhndl);
% %         clf;
% %     end
% %     set(fhndl,'Visible','on');  % assert visibility to bring window to front again after keypress
% %     drawFix(dat,dat.fix,[dat.I2MCopt.xres dat.I2MCopt.yres],img,[dat.I2MCopt.missingx dat.I2MCopt.missingy],sprintf('subj %s, trial %03d, stim: %s',files(p).subj,files(p).runnr,what{1}));
% %     pause
% %     if ~ishghandle(fhndl)
% %         return;
% %     end
% end
% if ishghandle(fhndl)
%     close(fhndl);
% end
% 
% rmpath(genpath(dirs.funclib));                  % cleanup path

