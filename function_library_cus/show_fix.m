function [startT, headDist] = show_fix(wpnt, x, y, fixTime, fixClrs, winRect, withTobii, tobiiFreq, EThndl, monWidth, monHeight)
    % To show a central fixation dot for fixTime or longer. 
    % If Tobbi is connected but the eyes are not tracked, 
    % or the fixation of eyes are not at the center, this 
    % function will wait 0.5s and check again, untill the 
    % eyes correctly fix to the center. 
    % 
    % Parameters
    % wpnt:      the window pointer
    % x:         x coor. of fixation center
    % y:         y coor. of fixation center
    % fixTime:   the minimum duration of fixation
    % fixClrs:   color of fixation
    % winRect:   window rect of Psychtoolbox
    % withTobii: skip steps in which tobii connection is necessary
    % tobiiFreq: sampling frequency of tobii 
    % EThndl:    the handle of a Titta object
    % monWidth:  physical width of monitor
    % monHeight: physical height of monitor
    %
    % Return
    % startT:    Psychtoolbox timestamp of fixation occurrence
    % headDist:  distance between head and screen, in centimeter 
    Screen('gluDisk',wpnt,fixClrs(1),x,y,round(winRect(3)/150));
    
    startT = Screen('Flip',wpnt);
    % log when fixation dot appeared in eye-tracker time. NB:
    % system_timestamp of the Tobii data uses the same clock as
    % PsychToolbox, so startT as returned by Screen('Flip') can be used
    % directly to segment eye tracking data 
    WaitSecs(fixTime); 
    while 1
        if withTobii
            SysbaseT = EThndl.getTimeAsSystemTime()-1e6;
            gazeData = EThndl.buffer.peekTimeRange('gaze',SysbaseT);
            headDists = getHeadDist(gazeData);
            headDist = headDists(3, find(~isnan(headDists(3, :)), 1, 'last'));
            if ~isempty(headDist)
                [startFixPix,~] = getLastFix(gazeData, monWidth, monHeight, tobiiFreq, headDist, winRect(3:4), SysbaseT, 15, 60);
                Fixerr = atan(norm(startFixPix-[winRect(3)/2, winRect(4)/2])/headDist);
                if Fixerr<1
                    break
                else
                    disp('Not fixed at central dot... try again...')
                end
            else 
                disp('EYE LOCATION is lost... try again...');
            end
        else
            headDist = 65; % cm, if head positioning was performed correctly. 
            break
        end
        checkend;
        WaitSecs(0.5);% reduce the sampling rate to reduce the pressure of CPU
        checkend;
    end
end