function [startT, headDist] = show_fix(wpnt, x, y, fixTime, fixClrs, winRect, MaxErr, varargin)
    % Show a fixation dot with support for both Tobii (TITTA) and EyeLink trackers
    
    % Parse input arguments
    p = inputParser;
    addParameter(p, 'eyeTrackerType', '', @ischar);
    addParameter(p, 'withTobii', false, @islogical);
    addParameter(p, 'tobiiFreq', 250, @isnumeric);
    addParameter(p, 'EThndl', [], @(x) isstruct(x) || isempty(x));
    addParameter(p, 'el', [], @(x) isstruct(x) || isempty(x));
    addParameter(p, 'monWidth', 51.1, @isnumeric);
    addParameter(p, 'monHeight', 28.7, @isnumeric);
    
    parse(p, varargin{:});
    
    eyeTrackerType = p.Results.eyeTrackerType;
    withTobii = p.Results.withTobii;
    tobiiFreq = p.Results.tobiiFreq;
    EThndl = p.Results.EThndl;
    el = p.Results.el;
    monWidth = p.Results.monWidth;
    monHeight = p.Results.monHeight;
    
    % Auto-detect eye tracker type
    if isempty(eyeTrackerType)
        if withTobii || (~isempty(EThndl) && isstruct(EThndl) && isfield(EThndl, 'buffer'))
            eyeTrackerType = 'Tobii';
        elseif ~isempty(el) && isstruct(el)
            eyeTrackerType = 'EyeLink';
        else
            eyeTrackerType = 'None';
        end
    end
    
    % Draw fixation dot 
    Screen('gluDisk', wpnt, fixClrs, x, y, round(winRect(3)/150));
    startT = Screen('Flip', wpnt);
    
    % Initialize
    headDist = 68;  % Default distance in cm
    fixationAcquired = false;
    acquireStartTime = GetSecs();
    
    % Device-specific fixation validation
    if strcmpi(eyeTrackerType, 'Tobii')
        %% TOBII: Wait and verify fixation at target location
        WaitSecs(fixTime);
        while ~fixationAcquired
            if ~isempty(EThndl) && isstruct(EThndl)
                SysbaseT = EThndl.getTimeAsSystemTime() - 1e6;
                gazeData = EThndl.buffer.peekTimeRange('gaze', SysbaseT);
                headDists = getHeadDist(gazeData);
                headDist = headDists(3, find(~isnan(headDists(3, :)), 1, 'last'));
                
                if ~isempty(headDist)
                    [startFixPix, ~] = getLastFix(gazeData, monWidth, monHeight, tobiiFreq, ...
                                                   headDist, winRect(3:4), SysbaseT, 15, 60);
                    % Calculate fixation error in degrees
                    Fixerr = atan(norm(startFixPix - [x, y]) / headDist) * 180 / pi;
                    
                    if Fixerr < MaxErr
                        fixationAcquired = true;
                    else
                        disp('Not fixed at central dot... try again...');
                    end
                else
                    disp('EYE LOCATION is lost... try again...');
                end
            else
                % No EThndl object available, just wait
                fixationAcquired = true;
            end
            
            if ~fixationAcquired
                checkend;
                WaitSecs(0.5);
                checkend;
            end
        end
        
    elseif strcmpi(eyeTrackerType, 'EyeLink')
        %% EYELINK: Wait and verify fixation at target location 
        % Continuously check gaze and wait until fixation is within MaxErr degrees
        
        fixationAcquired = false;
        fixStartTime = GetSecs();
        
        while ~fixationAcquired && (GetSecs() - fixStartTime) < fixTime
            checkend;
            
            % Read latest gaze sample
            if Eyelink('NewFloatSampleAvailable') > 0
                evt = Eyelink('NewestFloatSample');
                
                % Get which eye(s) are available
                eyeUsed = Eyelink('EyeAvailable');
                if eyeUsed == el.BINOCULAR
                    eyeUsed = el.LEFT_EYE;  % Default to left eye
                end
                
                % Extract gaze position in pixels based on available eye
                gazePos = [];
                if eyeUsed == el.LEFT_EYE && evt.gx(1) ~= -32768 && evt.gy(1) ~= -32768
                    gazePos = [evt.gx(1), evt.gy(1)];
                elseif eyeUsed == el.RIGHT_EYE && evt.gx(2) ~= -32768 && evt.gy(2) ~= -32768
                    gazePos = [evt.gx(2), evt.gy(2)];
                else
                    % Fallback: try any valid eye
                    if evt.gx(1) ~= -32768 && evt.gy(1) ~= -32768
                        gazePos = [evt.gx(1), evt.gy(1)];
                    elseif evt.gx(2) ~= -32768 && evt.gy(2) ~= -32768
                        gazePos = [evt.gx(2), evt.gy(2)];
                    end
                end
                
                % Validate fixation: check if gaze is within MaxErr degrees of target
                if ~isempty(gazePos)
                    % Convert pixel distance to degrees using same approach as Tobii
                    pixelDist = norm(gazePos - [x, y]);
                    % Angular error in degrees: arctan(pixel_distance / head_distance_in_pixels)
                    % Convert headDist from cm to pixels using monitor geometry
                    pixelsPerCm = winRect(3) / monWidth;  % pixels per cm
                    headDistPix = headDist * pixelsPerCm;
                    fixErr = atan(pixelDist / headDistPix) * 180 / pi;
                    
                    if fixErr < MaxErr
                        fixationAcquired = true;
                    else
                        % Not yet fixed at the target, continue checking
                    end
                end
            end
            
            WaitSecs(0.01);
        end
        
        if ~fixationAcquired
            disp('Not fixed at central dot... try again...');
        end
        
    else
        %% NO TRACKER: Just wait for the specified duration
        WaitSecs(fixTime);
        fixationAcquired = true;
    end
    
    if ~fixationAcquired
        disp('Warning: Fixation acquisition may have failed');
    end
end
