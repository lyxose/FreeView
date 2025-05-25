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
                                dirs.home       = cd;
cd data;                        dirs.data       = cd;
        cd samples_ophak;       dirs.samples    = cd;
cd ..;  cd fixDet;              dirs.fix        = cd;
cd ..;  cd msgs_ophak;          dirs.msgsO      = cd;
cd ..;  cd mat;                 dirs.mat        = cd;
cd ..;
cd ..;
cd function_library;            dirs.funclib    = cd;
cd ..;
cd results;                     dirs.res        = cd;
        cd 'AOImasks';          dirs.AOImasks   = cd;
cd ..;
cd(dirs.home);
addpath(genpath(dirs.funclib));                 % add dirs to path
addpath('function_library_cus');

binSize = 25; % in pixel
show_nFix = 8; % only plot the first n fixations

%%% get result table with eye data
[files,~] = FileFromFolder(dirs.fix,[],'mat');
filtstr = '^Dat_Sub(\d+)_Ses(\d+).mat$';
matched = regexpi({files.name}.',filtstr,'tokens');
resfile_idx = ~cellfun(@isempty,matched);
resfiles   = files(resfile_idx);
sub_ses_res = cellfun(@(x) [str2double(x{1}{1}), str2double(x{1}{2})], matched(resfile_idx), 'UniformOutput', false);
sub_ses_res = cell2mat(sub_ses_res(~cellfun(@isempty, matched(resfile_idx))));

select_sess = 1:length(resfiles);
% select_sess = 10:14; % single subj

alpha = 0.02; % traj

figure(1);  % 在循环外创建 figure，固定窗口
figure(2);  % 在循环外创建 figure，固定窗口

if ~exist('./results/heatmap/','dir')
    mkdir('./results/heatmap/');
end
if ~exist('./results/trajec/','dir')
    mkdir('./results/trajec/');
end


for k = 24:-1:1
    fixPos = [];
    stTs = [];   % color will be depended on the start time 
    triSpl = []; % to split each trial by the end index
    
for p=select_sess
    sess = load(fullfile(dirs.mat,sprintf("Dat_Sub%.0f_Ses%.0f.mat",sub_ses_res(p,1),sub_ses_res(p,2))),'expt','geometry');
    resT = load(fullfile(dirs.fix,[resfiles(p).fname '.mat'])).expT;

    [c,~,ic] = unique([resT.ECC, resT.Orient],'rows');
    img_width = sess.expt.winRect(3);  % 二维空间的宽度
    img_height = sess.expt.winRect(4); % 二维空间的高度
    edges_x = linspace(1, img_width,  round(img_width/binSize));  % 网格的X边界
    edges_y = linspace(1, img_height, round(img_height/binSize)); % 网格的Y边界
    ut = UT(sess.geometry.displayArea.width/10, img_width, mean(resT.headDist));
    tgWidth = ut.deg2pix(sess.expt.GaborWidth);
    target_loc = ut.Pol2Rect(c(k,:));
    target_loc = target_loc.*[1,-1]+[img_width, img_height]/2;
%     for k = 1:max(ic)
%         fixPos = [];
        iresT = resT(ic == k,:); % subset of resT, all in the EO ECC-Ori condition 
        for i=1:height(iresT)
            tFixPos = transpose([iresT.dat(i).fix.xpos; iresT.dat(i).fix.ypos]); % this trial eye trajectory
%             disp(norm(tFixPos(1,:)-[img_width, img_height]/2));
            tTime = iresT.dat(i).fix.startT;  % start time stamp
            tFixPos = tFixPos(tTime>0,:);
            tTime = tTime(tTime>0);
            while length(tFixPos(:,1))>1 && norm(tFixPos(1,:)-[img_width, img_height]/2)<2*binSize
                if length(tFixPos(:,1))==1
                    break
                end
                tFixPos = tFixPos(2:end,:); % remove the start fixation
                tTime = tTime(2:end,:);
            end
            % cut to early fixations
            if length(tTime(:,1))>show_nFix
                tFixPos = tFixPos(1:show_nFix,:);
                tTime = tTime(1:show_nFix,:);
            end
%             tFixPos = tFixPos(1,:);
            tsp(1)  = length(stTs)+1;
            stTs = [stTs;tTime];
            fixPos = [fixPos; tFixPos];
            tsp(2) = length(stTs);
            triSpl = [triSpl;tsp];
            
        end
end


        figure(1);
        clf;  
        counts = hist3(fixPos, 'Edges', {edges_x, edges_y});
        density_map = counts' / sum(counts(:));
%         figure;
        imagesc(edges_x, edges_y, density_map);  % 使用imagesc绘制热图
        pr = prctile(density_map(:), 99.8);
        clim([min(density_map(:)), pr]);
%         set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        axis image;
        colorbar; % 显示颜色条
        % 标记中央位置：用 rectangle 函数绘制一个圆形
        hold on;  % 保持当前图形
        rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
                  'EdgeColor', '#800080', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
        rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
                  'EdgeColor', 'r', 'LineWidth', 0.7, 'Curvature', [1, 1]);  %        

        title(sprintf('Fixation Heatmap %.0f-ECC %.0f-Ori  %.0f subj',c(k,1),c(k,2),length(select_sess)));
%         set(gcf, 'PaperPositionMode', 'auto');
        exportgraphics(gcf, sprintf('./results/heatmap/%.0fECC_%.0fOri_%.0fsess.png',c(k,1),c(k,2),length(select_sess)), 'Resolution', 600);
        savefig(sprintf('./results/heatmap/%.0fECC_%.0fOri_%.0fsess.fig',c(k,1),c(k,2),length(select_sess)))
%         print(gcf,sprintf('./results/heatmap/%.0fECC_%.0fOri_%.0fsubj.png',c(k,1),c(k,2),length(select_sess)), '-dpng', '-r600', '-tight')
        hold off;

        figure(2); clf; hold on; axis equal; 
        colormap(cool);  % 选择traj的颜色映射，例如 'parula', 'jet', 'hot'
        clim([0, 1]);
    
        % 归一化时间值到 [0,1]，用于颜色映射
%             normTime = time / maxTime;
%         normTime = tTime / max(tTime);
        rankTime = tiedrank(stTs);
        normTime = rankTime/max(rankTime);
        
        % 绘制渐变线段
        for i = 1:size(triSpl,1)
            tsp = triSpl(i,:);
            xpos = fixPos(tsp(1):tsp(2),1);
            ypos = fixPos(tsp(1):tsp(2),2);
            tTime = normTime(tsp(1):tsp(2));
            for j = 1:length(xpos)-1
                ct = tTime(j);  % 当前时间的归一化值
                color = cool(256);  % 获取 colormap 颜色
                idx = max(1, round(ct * 255));  % 找到颜色索引
                plot(xpos(j:j+1), ypos(j:j+1), '-', 'Color', [color(idx, :), alpha], 'LineWidth', 1.5);
            end
        end
        % 设置图像的坐标范围和 y 轴的零点在左上角
        axis([0 img_width 0 img_height]);  % 横纵坐标范围为 img_width 和 img_height
        set(gca, 'YDir', 'reverse');  % y轴零点在左上角
        % 绘制目标圆圈
        rectangle('Position', [round(img_width/2)-5, round(img_height/2)-5, 10, 10], ...
                  'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);  % 红色圆形
        rectangle('Position', [target_loc(1)-tgWidth/2, target_loc(2)-tgWidth/2, tgWidth, tgWidth], ...
          'EdgeColor', 'k', 'LineWidth', 0.7, 'Curvature', [1, 1]);
        colorbar;  % 显示颜色条
        title(sprintf('Fixation Trajectory %.0f-ECC %.0f-Ori  %.0f subj',c(k,1),c(k,2),length(select_sess)));
%         xlabel('X Position');
%         ylabel('Y Position');
        exportgraphics(gcf, sprintf('./results/trajec/%.0fECC_%.0fOri_%.0fsess_tra.png',c(k,1),c(k,2),length(select_sess)), 'Resolution', 600);
        savefig(sprintf('./results/trajec/%.0fECC_%.0fOri_%.0fsess_tra.fig',c(k,1),c(k,2),length(select_sess)))
        hold off; 
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
