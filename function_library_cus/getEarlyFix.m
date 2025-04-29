%% ophakker    

function [earlyFixPix, dat] = getEarlyFix(gazeData, firstn, monWidth, monHeight, tobiiFreq, headDist, scrRes, SysstartT, maxMergeDist, minFixDur)
    
    if 0
        firstn=2; % first 2 fixation will return
        monWidth=54.47853; % EThndl.geom.displayArea.width, in mm
        monHeight=30.64417; %EThndl.geom.displayArea.height, in mm
        tobiiFreq = 250;
        headDist = 68;
        SysstartT = 1282425144;
        scrRes  = [3480,2160];% dat.expt.resolution;
        maxMergeDist = 15;
        minFixDur = 60;
    end
    
    % the Pro SDK does not guarantee invalid data is nan. Set to nan if
    % invalid
    ts      = gazeData.systemTimeStamp;
    gazeData. left.gazePoint.onDisplayArea(:,~gazeData. left.gazePoint.valid) = nan;
    gazeData.right.gazePoint.onDisplayArea(:,~gazeData.right.gazePoint.valid) = nan;
    gazeData. left.pupil.diameter(~gazeData. left.pupil.valid) = nan;
    gazeData.right.pupil.diameter(~gazeData.right.pupil.valid) = nan;
    % collect data from the file, and turn gaze positions from normalized
    % coordinates into pixels
    samp    = [bsxfun(@times,gazeData.left.gazePoint.onDisplayArea,scrRes.'); bsxfun(@times,gazeData.right.gazePoint.onDisplayArea,scrRes.'); gazeData.left.pupil.diameter; gazeData.right.pupil.diameter];
        
        
    % data
    data = transpose(double([ts; samp]));
    
    
    %% I2MC
        
    % load messges and trial mat file. We'll need to find when in trial the
    % stimulus came on to use that as t==0
    
    % event detection
    % make params struct (only have to specify those you want to be
    % different from their defaults)
    opt.xres          = scrRes(1);
    opt.yres          = scrRes(2);
    opt.missingx      = nan;
    opt.missingy      = nan;
    opt.scrSz         = [monWidth monHeight];  % cm, EThndl.geom.displayArea.width
    opt.disttoscreen  = headDist;
    opt.freq          = tobiiFreq;
    if opt.freq>120
        opt.downsamples   = [2 5 10];
        opt.chebyOrder    = 8;
    elseif opt.freq==120
        opt.downsamples   = [2 3 5];
        opt.chebyOrder    = 7;
    else
        % 90 Hz, 60 Hz, 30 Hz
        opt.downsampFilter= false;
        opt.downsamples   = [2 3];
    end
    if opt.freq==30
        warning('Be careful about using I2MC with data that is only 30 Hz. In a brief test, this did not appear to work well with the settings in this file.')
    end
    if (~isfield(opt,'downsampFilter') || opt.downsampFilter) && ~exist('cheby1','file')
        warning('By default, I2MC runs a Chebyshev filter over the data as part of its operation. It appears that this filter (the function ''cheby1'' from the signal processing toolbox) is not available in your installation. I am thus disabling the filter.')
        opt.downsampFilter= false;
    end
    opt.maxMergeDist  = maxMergeDist;
    opt.minFixDur     = minFixDur;
    
    % make data struct
    clear dat;
    dat.time        = (data(:,1)-double(SysstartT))./1000; % mu_s to ms, make samples relative to onset of picture
    dat.left.X      = data(:,2);
    dat.left.Y      = data(:,3);
    dat.right.X     = data(:,4);
    dat.right.Y     = data(:,5);
    dat.left.pupil  = data(:,6);    % add pupil data to file. not used by I2MC but good for plotting
    dat.right.pupil = data(:,7);
    [fix,dat]       = I2MCfunc(dat,opt);
    
    % collect info and store
    dat.fix         = fix;
    dat.I2MCopt     = opt;
    %% last fixation location
    earlyFixPix = [transpose(dat.fix.xpos(1:firstn)), transpose(dat.fix.ypos(1:firstn))]; % [0,0] for upperleft corner
    