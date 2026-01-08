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
    Screen('gluDisk', wpnt, fixClrs(1), x, y, round(winRect(3)/150));
    startT = Screen('Flip', wpnt);
    
    % Initialize
    headDist = 68;  % Default distance
    fixationAcquired = false;
    startTime = GetSecs();
    
    % Common fixation validation loop for all tracker types
    while ~fixationAcquired   % no timeout
        
        % Get gaze data in unified format: [x, y, headDist]
        gazePos = [];
        
        if strcmpi(eyeTrackerType, 'Tobii')
            if ~isempty(EThndl) && isstruct(EThndl)
                SysbaseT = EThndl.getTimeAsSystemTime() - 1e6;
                gazeData = EThndl.buffer.peekTimeRange('gaze', SysbaseT);
                headDists = getHeadDist(gazeData);
                headDist = headDists(3, find(~isnan(headDists(3, :)), 1, 'last'));
                
                if ~isempty(headDist)
                    [startFixPix, ~] = getLastFix(gazeData, monWidth, monHeight, tobiiFreq, ...
                                                   headDist, winRect(3:4), SysbaseT, 15, 60);
                    gazePos = startFixPix;  % [x, y]
                else
                    disp('EYE LOCATION is lost... try again...');
                end
            end
            
        elseif strcmpi(eyeTrackerType, 'EyeLink')
            if ~isempty(el)
                % Get latest eye data from EyeLink
                evt = Eyelink('GetNextDataFile');
                if evt.type == el.FIXATION
                    gazePos = [evt.gavx, evt.gavy];  % Gaze position in pixels
                    % Optionally get headDist from EyeLink if available
                end
            end
        end
        
        % Common validation logic for both trackers
        if ~isempty(gazePos)
            % Calculate fixation error (angular deviation in radians)
            Fixerr = atand(norm(gazePos - [x, y]) / headDist);
            
            % Check if fixation is on target (< 1 degree threshold)
            if Fixerr < MaxErr
                % Verify fixation duration
                if (GetSecs() - startTime) >= fixTime
                    fixationAcquired = true;
                end
            else
                disp('Not fixed at central dot... try again...');
            end
        end
        
        if ~fixationAcquired
            checkend;
            WaitSecs(0.1);  % Brief pause before next check
            checkend;
        end
    end
    
    if ~fixationAcquired
        disp('Warning: Fixation not acquired within timeout period');
    end
end
