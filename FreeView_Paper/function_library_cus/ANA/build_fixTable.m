function [fixTable, Nsubj, img_width, img_height, ut, center, digPlace] = build_fixTable(select_sess, exclude_sub, sub_ses_res, dirs, resfiles, learn_stage_n, last_trial, skip_corr)
    % 构建注视点数据表格
    % 输入：
    % select_sess: 要处理的会话索引列表
    % exclude_sub: 要排除的被试ID列表
    % sub_ses_res: 被试-会话对应表
    % dirs: 目录结构体，包含.mat和fix目录路径
    % resfiles: 结果文件信息结构体数组
    % learn_stage_n: 学习阶段的试次数
    % last_trial: 最后一个试次编号
    % skip_corr: 是否跳过正确响应的试次（布尔值）
    % 输出：
    % fixTable: 注视点数据表格
    % Nsubj: 处理的被试数量
    % img_width, img_height: 图像宽高
    % ut: 坐标转换工具对象
    % center: 屏幕中心坐标
    % digPlace: 数字位置数据（如果有）
    fixTable = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[], ...
        'VariableNames',{'subID','sessID','Ntrial','xpos','ypos','r','theta','tgX','tgY','tgr','tgtheta','startT','dur','Nfix','lNfix'});
    Nsubj=0;
    img_width = [];
    img_height = [];
    ut = [];
    center = [];
    for p=select_sess
        if ismember(sub_ses_res(p,1),exclude_sub)
            continue
        end
        Nsubj=Nsubj+1;
        sess = load(fullfile(dirs.mat,sprintf("Dat_Sub%.0f_Ses%.0f.mat",sub_ses_res(p,1),sub_ses_res(p,2))),'expt','geometry');
        resT = load(fullfile(dirs.fix,[resfiles(p).fname '.mat'])).expT;
        try
            img_width = sess.expt.winRect(3);
            img_height = sess.expt.winRect(4);
        catch
            img_width=1920;
            img_height=1080;
        end
        center = [img_width,img_height]./2;
        ut = UT(sess.geometry.displayArea.width/10, img_width, mean(resT.headDist,1, 'omitnan'), false);
        iresT = resT;
        % 如果存在 dirs.thrData，则读取阈值阶段的表格数据并拼接
        if isfield(dirs, 'thrDir') && exist(dirs.thrDir, 'dir')
            if sub_ses_res(p,1) == 17
                % 17号被试需要读取 ses0 和 ses1，0在前，1在后
                thrT0 = load(fullfile(dirs.thrDir, sprintf("Dat_Sub17_Ses0.mat")));
                thrT1 = load(fullfile(dirs.thrDir, sprintf("Dat_Sub17_Ses1.mat")));
                thrTab = [thrT0.expT; thrT1.expT];
            else
                thrTab = load(fullfile(dirs.thrDir, sprintf("Dat_Sub%.0f_Ses%.0f.mat", sub_ses_res(p,1), sub_ses_res(p,2)))).expT;
            end
            % 拼接到 iresT 前面（纵向），字段名需一致
            iresT = [thrTab; iresT];
        end
        for i=1:height(iresT)
            if i <= learn_stage_n || i>last_trial
                continue
            end
            if skip_corr && ~isnan(iresT.key2RT(i)) || isempty(iresT.dat(i).fix)
                continue
            end
            idx = iresT.dat(i).fix.startT(:)>0;
            tX = iresT.dat(i).fix.xpos(idx)';
            tY = iresT.dat(i).fix.ypos(idx)';
            tStartT = iresT.dat(i).fix.startT(idx);
            tDur = iresT.dat(i).fix.dur(idx);
            nFixs = numel(tX);
            if nFixs>0
                ttgLoc = ut.deg2pix([iresT.Xtarg(i),iresT.Ytarg(i)]) .* [1,-1] + center;
                ttgX = ttgLoc(1);
                ttgY = ttgLoc(2);
                ttgR = iresT.ECC(i);
                ttgTheta = iresT.Orient(i);
                subID = repmat(sub_ses_res(p,1), nFixs, 1);
                sessID = repmat(sub_ses_res(p,2), nFixs, 1);
                Ntrial = repmat(i, nFixs, 1);
                Nfix = (1:nFixs)';
                lNfix = (nFixs:-1:1)';
                trial_rTheta=ut.Rect2Pol(([tX,tY]-center).*[1,-1]);
                tgX = repmat(ttgX, nFixs, 1);
                tgY = repmat(ttgY, nFixs, 1);
                tgr = repmat(ttgR, nFixs, 1);
                tgtheta = repmat(ttgTheta, nFixs, 1);
                newRows = table(subID, sessID, Ntrial, tX, tY, trial_rTheta(:,1), trial_rTheta(:,2),tgX, tgY, tgr, tgtheta, tStartT, tDur, Nfix, lNfix,  ...
                    'VariableNames',{'subID','sessID','Ntrial','xpos','ypos','r','theta','tgX','tgY','tgr','tgtheta','startT','dur','Nfix','lNfix'});
                fixTable = [fixTable; newRows];
            end
        end
    end
    fixTable.Err = sqrt((fixTable.xpos - fixTable.tgX).^2 + (fixTable.ypos - fixTable.tgY).^2);
    fixTable.degErr = ut.pix2deg(fixTable.Err);
    saccades = [0,0;ut.Rect2Pol(([diff(fixTable.xpos), diff(fixTable.ypos)]).* [1,-1])];
    saccades(fixTable.Nfix==1,:) = 0;
    fixTable.saccLen   = saccades(:,1);
    fixTable.saccAngle = saccades(:,2);
    digPlace = floor(log10(max(fixTable.Ntrial)));
    fixTable.TriID = fixTable.subID*10^(digPlace+2) + fixTable.sessID*10^(digPlace+1) + fixTable.Ntrial;
end
