function saveStiImg(item, seed, path)
    if nargin<2
        seed = randi(1000);
    end
    if nargin<3
        path='';
    end
    
    winRect(3) = 1920;
    winRect(4) = 1080;
    ut = UT(40, winRect(3), 65);
    if strcmp(item,'Gabor')
        bgContrast = 0;
        tgContrast = 1;
        GaborWidth = 10;
        GaborSF = 2/GaborWidth;
        GaborOrient = -45;
        bgWidth = 15;
        tgCenter = [0,0];
    elseif strcmp(item,'bg')
        bgContrast = 0.2;
        tgContrast = 0;
        GaborWidth = 0;
        GaborSF = 2/GaborWidth;
        GaborOrient = -45;
        bgWidth = 15;
        tgCenter = [0,0];
    elseif strcmp(item,'T1')
        bgContrast = 0.15;
        tgContrast = 0.2;
        GaborWidth = 1/3;
        GaborSF = 2/GaborWidth;
        GaborOrient = -45;
        bgWidth = 15;
        tgCenter = [0,0];
    elseif strcmp(item,'T1evid')
        bgContrast = 0.15;
        tgContrast = 0.3;
        GaborWidth = 2/3;
        GaborSF = 2/GaborWidth;
        GaborOrient = -45;
        bgWidth = 15;
        tgCenter = [0,0];
    elseif strcmp(item,'T2')
        bgContrast = 0.15;
        tgContrast = 0.2;
        GaborWidth = 1/3;
        GaborSF = 2/GaborWidth;
        GaborOrient = -45;
        bgWidth = 15;
        tgCenter = ut.Pol2Rect([4,45]);
    elseif strcmp(item,'T2evid')
        bgContrast = 0.15;
        tgContrast = 0.3;
        GaborWidth = 2/3;
        GaborSF = 2/GaborWidth;
        GaborOrient = -45;
        bgWidth = 15;
        tgCenter = ut.Pol2Rect([4,45]);
    elseif strcmp(item,'T2space') % with 24 possible position
        bgContrast = 0.15;
        tgContrast = 0.3;
        GaborWidth = 2/3;
        GaborSF = 2/GaborWidth;
        GaborOrient = -45;
        bgWidth = 15;
        tgCenter = ut.Pol2Rect([4,45]);
    end
    stimulus = genStim(winRect, ut, bgContrast, tgContrast, tgCenter, GaborSF, GaborWidth, GaborOrient, bgWidth, seed);
    figure();
    imshow(stimulus);
    hold on
    if strcmp(item,'T2space')
        % 标记目标位置
        for ecc = [2,4,6]
        for ori = 0:45:315
            tgWidth = ut.deg2pix(GaborWidth);
            target_loc = ut.Pol2Rect([ecc,ori]);
            target_loc = target_loc.*[1,-1]+winRect(3:4)/2;
            rectangle('Position', [target_loc(1)-tgWidth*0.8, target_loc(2)-tgWidth*0.8, 1.6*tgWidth, 1.6*tgWidth], ...
                      'EdgeColor', '#E3170D', 'LineWidth', 0.6, 'Curvature', [1, 1]);
        end
        end
    end
    hold off
    exportgraphics(gcf, [path,'/',item,'.png'], 'Resolution', 300);  % 300 DPI保存图像