% Functions in Experiment Run
classdef expRun
    methods (Static)
        %%%%%%%%%%%%%%%%%%%% Preparation Functions %%%%%%%%%%%%%%%%%%%%%%
        function obj = setProperties(obj, varargin)
            % set multiple properties using varargin, which is a cell array of input arguments
            propArgs = varargin;
            for i = 1:2:length(propArgs)
                propName = propArgs{i};
%                 disp(propName)
                propValue = propArgs{i+1};
                obj.(char(propName)) = propValue;
            end
        end 
        %------------------------------
        function [trialList,trialTable] = generateTrialList(varargin)
        % generate trial list
        % orthogonal design
            if rem(nargin,2)~=0
                error('Incorrect Input Format! Correct should be: VariableName,Variable(x*1),...')
            else
                fprintf('Num of Input Variables: %d\n',nargin)
                % celldisp(varargin)
                nvar = nargin/2;
                varNames = cell(1,nvar);
                varValue = cell(1,nvar);
                dim = 1;
                for i = 1:nvar
                    varNames{i} = varargin{i*2-1};
                    varValue{i} = varargin{i*2};
                    if ~iscell(varValue{i})
                        varValue{i} = num2cell(varValue{i});
                    end
                    if ~iscolumn(varValue{i})
                        if isrow(varValue{i})
                            varValue{i} = varValue{i}';
%                             disp(varValue{i})
                        else
                            error('Variable should be organized in a vector!!') 
                        end
                    end
                    dim = dim*length(varValue{i});
                end
                % create trial list as a table
                trialList = cell(dim,nvar);
                n1 = length(varValue{1});
                trialList(:,1) = repmat(varValue{1},dim/n1,1);
                for i = 2:nvar
                    n2 = n1*length(varValue{i});
                    trialList(:,i) = repmat(repelem(varValue{i},n1,1),dim/n2,1);
                    n1 = n1*length(varValue{i});
                end        
                trialTable = cell2table(trialList);
                trialTable.Properties.VariableNames = varNames;         
            end
        end
        %--------------------------------
        function key = keyChose(funkey,respkey)
            % Prepare the keys used
            KbName('UnifyKeyNames');
            key.esc         = KbName(funkey{1}); 
            key.continue    = KbName(funkey{2});
            key.alt         = cell(1,length(respkey));
            for i = 1:length(respkey)
                key.alt{i}  = KbName(respkey{i});
            end
        end
        
        %%%%%%%%%%%%%%%%%%%% Screen Setup Functions %%%%%%%%%%%%%%%%%%%
        function [win,rect,ifi,center] = openScreen
            % prepare windows
            HideCursor;
            AssertOpenGL;
            whichscreen = max(Screen('Screens'));

            black = BlackIndex(whichscreen);
            white = WhiteIndex(whichscreen);
            gray = round((black+white)/2);

            [win, rect] = Screen('OpenWindow',whichscreen,gray);
            [center(1), center(2)] = RectCenter(rect);

            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            
            ifi = Screen('GetFlipInterval',win);

            % open window
            Screen('TextSize',win, 30);
            DrawFormattedText(win,'Welcome to our experiment! \n \n Press the continue key to start!'...
                    ,'center','center',200);
            Screen('Flip', win);
        end
        %---------------------------
        function closeScreen
             ShowCursor; 
             Screen('CloseAll');   
        end
        
        %---------------------------
        function PPM = cardTask(funkey,respkey)
            % specify key
            % respkey towards (left,right,up,down)
            key = expRun.keyChose(funkey,respkey);
            % open screen
            [win, ~, ~ ,center] = expRun.openScreen;
            % begin test 
            expRun.keyPressToContinue(key);
            % present virtual card
            % initialize card rect and stepsize
            cWidth      = 300;
            cHeight     = 200;
            cardrect    = [center(1)-cWidth/2,center(2)-cHeight/2,center(1)+cWidth/2,center(2)+cHeight/2]; 
            stp         = 1; % in pixel
            color       = 0;
            pressDelay  = 0.05; % in secs
            % key press check
            [~, ~, keycode] = KbCheck;
            while ~keycode(key.esc)
                [~, ~, keycode] = KbCheck;
                Screen('FillRect', win, color, cardrect);
                Screen('Flip',win);
                if keycode(key.alt{1})==1
                    cardrect(3) = cardrect(3)-stp;
                elseif keycode(key.alt{2})==1
                    cardrect(3) = cardrect(3)+stp;
                elseif keycode(key.alt{3})==1
                    cardrect(4) = cardrect(4)-stp;
                elseif keycode(key.alt{4})==1
                    cardrect(4) = cardrect(4)+stp;
                end
                WaitSecs(pressDelay);
            end
            % card physical distance 85 mm (default)
            % calculate pixel per mm (PPM)
            PPM = (cardrect(3)-cardrect(1))/85;  
            % close screen
            expRun.closeScreen;
        end
        %-------------------------------
        function [viewDistance, dotXY] = blindSpotTask(funkey,respkey,PPM)
            % PPM calculated from cardtask or experimenter 
            % specify key
            % respkey towards (left,right)
            key = expRun.keyChose(funkey,respkey);
            % open screen
            [win, ~, ~, center] = expRun.openScreen;
            % begin test 
            expRun.keyPressToContinue(key);
            % initialize params
            dotDia      = 20; % in pixel
            dotSpeed    = 2; % in pixel
            dotXY       = [-300,0];
            dotColor    = [255 0 0]; % RGB 
            pressDelay  = 0.05; % in secs
            fixLen      = 10;
            fixColor    = 0;
            % key press check
            [~, ~, keycode] = KbCheck;
            while ~keycode(key.esc)
                [~, ~, keycode] = KbCheck;
                Screen('DrawDots', win, dotXY, dotDia, dotColor, center,1);
                expRun.crossFixation(win,center,fixLen,fixColor,2);
                Screen('Flip',win);
                if keycode(key.alt{1})==1
                    dotXY(1) = dotXY(1)-dotSpeed;
                elseif keycode(key.alt{2})==1
                    dotXY(1) = dotXY(1)+dotSpeed;
                end
                WaitSecs(pressDelay);
            end
            % visual angle of blind spot: 13.5 degree (default, according to Li et al., 2020, sci reports)
            % calculate view distance
            viewDistance = abs(dotXY(1))/PPM/tan(13.5/180*pi);  
            % close screen
            expRun.closeScreen;
            
        end        
        %--------------------------       
        function [PPM,viewDistance] = setupScreen
            % rect: screen resolution (e.g., [1024,768]), only the width is
            % used in the calculation below
            %--------------
            % dlg input
            dlgprompt = {'Sreen Width (mm):',...
                'View Distance (mm):'};
            dlgname       = 'expSetup';
            numlines      = 1;
            defaultans    ={'',''};
            expSetup      = inputdlg(dlgprompt,dlgname,numlines,defaultans);
            %---------------
            % pixel per mm
            rect = get(0,'screensize');
            if ~isempty(expSetup{1})
                PPM = rect(3)./str2double(expSetup{1});
            else
                funKey = {'escape','space'};
                respKey = {'LeftArrow','RightArrow','UpArrow','DownArrow'};
                fprintf('Card Task Begin !! \n');
                PPM = expRun.cardTask(funKey,respKey);
            end
            fprintf('PPM: %0.2f \n',PPM);
            WaitSecs(2);
                    
            % view distance
            if ~isempty(expSetup{2})
                viewDistance = str2double(expSetup{2});
            else
                funKey = {'escape','space'};
                respKey = {'LeftArrow','RightArrow'};
                fprintf('Blind Spot Task Begin !! \n');
                viewDistance = expRun.blindSpotTask(funKey,respKey,PPM);
            end
            fprintf('View Distance: %0.2f \n',viewDistance);
            WaitSecs(2);
        end
        %-----------------------
        function P = angle2pixel(A,PPM,viewDistance)
            % from angle to pixel
            P = tan(A).*viewDistance.*PPM;            
        end
        %------------------------
        function A = pixel2angle(P,PPM,viewDistance)
            % from pixel to angle
            A = atan(P./viewDistance./PPM);            
        end
        
        %%%%%%%%%%%%%%%%%%%% Experiment Run Functions %%%%%%%%%%%%%%%
        function dispInstruction(win,expName)
            % display instruction   
            insImg = read(strcat('Instruction/',expName));
            insTex = Screen('MakeTexture', win, insImg);
            Screen('DrawTexture', win,insTex); 
            Screen('Flip', win);
        end
        %-------------------------
        function [toContinue, toQuit] = keyPressToContinue(key)
            % waiting for Key Pressed To Continue/Quit
            [~, ~, keycode] = KbCheck;
            while keycode(key.continue)==0 && keycode(key.esc)==0
              [~, ~, keycode] = KbCheck;
            end
            toContinue = keycode(key.continue);
            toQuit = keycode(key.esc);
        end
        %----------------------------
        function crossFixation(win,center,fixLen,fixColor,penWidth)
            % fixation
            fixXY = [0 0 -fixLen fixLen;-fixLen fixLen 0 0]; % a cross
            Screen('Drawlines', win, fixXY, penWidth, fixColor, center);              
        end
        %-----------------------------
        function drawPeriCues(win,cueDegree,cueEccentricity,targDia,cueColor,center,penWidth,fillOval) 
            % cueDegree, cueEccentricity, targDia, penWidth, and center in pixels
            % cueDegree should be a column vector, multiple cues would be drawn then.
            % cueEccentricity should be a scalor or a column vector with the size of cueDegree.
            % transform the polar cue xy to cartesian         
            nCue = length(cueDegree);
            if isscalar(cueEccentricity)
                cueEccentricity = ones(nCue,1).*cueEccentricity;             
            end            
            [cueCenter(1,:), cueCenter(2,:)] = pol2cart(cueDegree,cueEccentricity);
            cueRect = [cueCenter(1,:)-targDia/2;cueCenter(2,:)-targDia/2;cueCenter(1,:)+targDia/2;cueCenter(2,:)+targDia/2];
            if fillOval==0
                Screen('FrameOval', win, cueColor, cueRect+repmat(center',2,nCue), penWidth);
            else
                Screen('FillOval', win, cueColor, cueRect+repmat(center',2,nCue), penWidth);
            end
        end  
        %------------------------------
        function drawArrow(win, cueDeg, arrowLen, cueColor, center, penWidth, distanceFromCenter)
            % arrowLen: the length of the arrow wings
            % cueDeg: cue angle (in radians)
            % cueColor: color of the arrow
            % center: coordinates of the screen center
            % penWidth: line width
            % distanceFromCenter: distance between the arrow wings and the fixation point

            % Define the angle offset for the wings
            wingAngle = pi / 6; % 30 degrees offset

            % Preliminary wing point positions based on distance from center
            tempWing1X = center(1) + distanceFromCenter * cos(cueDeg);
            tempWing1Y = center(2) - distanceFromCenter * sin(cueDeg);
            tempWing2X = center(1) + distanceFromCenter * cos(cueDeg + pi); % opposite direction
            tempWing2Y = center(2) - distanceFromCenter * sin(cueDeg + pi); % opposite direction

            % Calculate left and right endpoints for each wing
            % First wing
            wing1LeftX = tempWing1X - arrowLen * cos(cueDeg - wingAngle);
            wing1LeftY = tempWing1Y + arrowLen * sin(cueDeg - wingAngle);
            wing1RightX = tempWing1X - arrowLen * cos(cueDeg + wingAngle);
            wing1RightY = tempWing1Y + arrowLen * sin(cueDeg + wingAngle);

            % Second wing
            wing2LeftX = tempWing2X - arrowLen * cos(cueDeg - wingAngle);
            wing2LeftY = tempWing2Y + arrowLen * sin(cueDeg - wingAngle);
            wing2RightX = tempWing2X - arrowLen * cos(cueDeg + wingAngle);
            wing2RightY = tempWing2Y + arrowLen * sin(cueDeg + wingAngle);

            % Calculate the midpoint of each wing segment
            wing1MidX = (wing1LeftX + wing1RightX) / 2;
            wing1MidY = (wing1LeftY + wing1RightY) / 2;
            wing2MidX = (wing2LeftX + wing2RightX) / 2;
            wing2MidY = (wing2LeftY + wing2RightY) / 2;

            arrow1MidX = (wing1MidX + tempWing1X)/2;
            arrow1MidY = (wing1MidY + tempWing1Y)/2;    
            arrow2MidX = (wing2MidX + tempWing2X)/2;
            arrow2MidY = (wing2MidY + tempWing2Y)/2;

            % Calculate the midpoint between the two wing midpoints
            arrowCenterX = (arrow1MidX + arrow2MidX) / 2;
            arrowCenterY = (arrow1MidY + arrow2MidY) / 2;

            % Adjust the calculated arrow center to the desired distance from center
            offsetX = center(1) - arrowCenterX;
            offsetY = center(2) - arrowCenterY;
            Wing1X = tempWing1X + offsetX;
            Wing1Y = tempWing1Y + offsetY;
            Wing2X = tempWing2X + offsetX;
            Wing2Y = tempWing2Y + offsetY;

            % Adjusted positions for first arrow wing with offset
            wing1LeftX = Wing1X - arrowLen * cos(cueDeg - wingAngle);
            wing1LeftY = Wing1Y + arrowLen * sin(cueDeg - wingAngle);
            wing1RightX = Wing1X - arrowLen * cos(cueDeg + wingAngle);
            wing1RightY = Wing1Y + arrowLen * sin(cueDeg + wingAngle);

            % Draw first wing
            Screen('DrawLine', win, cueColor, Wing1X, Wing1Y, wing1LeftX, wing1LeftY, penWidth);
            Screen('DrawLine', win, cueColor, Wing1X, Wing1Y, wing1RightX, wing1RightY, penWidth);

            % Adjusted positions for second arrow wing with offset
            wing2LeftX = Wing2X - arrowLen * cos(cueDeg - wingAngle);
            wing2LeftY = Wing2Y + arrowLen * sin(cueDeg - wingAngle);
            wing2RightX = Wing2X - arrowLen * cos(cueDeg + wingAngle);
            wing2RightY = Wing2Y + arrowLen * sin(cueDeg + wingAngle);

            % Draw second wing
            Screen('DrawLine', win, cueColor, Wing2X, Wing2Y, wing2LeftX, wing2LeftY, penWidth);
            Screen('DrawLine', win, cueColor, Wing2X, Wing2Y, wing2RightX, wing2RightY, penWidth);
        end
        %------------------------------
        function [lineL, lineT, lineFrame] = LTtask
            % generate the line coordinates for L, T and frame
            % center (0,0)
            % lines on the frame (3 horz 3 vert)
            [x,y] = meshgrid(-1:1,-1:1);             
            x = reshape([x(:,[1,3])',x([1,3],:)],1,2*6);
            y = reshape([y(:,[1,3])',y([1,3],:)],1,2*6);
            lineFrame = [x;y];
            % rotate
            lineL(:,:,1) = reshape(lineFrame(:,[1,2,7,8]),2,2*2); % 0
            lineL(:,:,2) = reshape(lineFrame(:,[1,2,11,12]),2,2*2); % 90
            lineL(:,:,3) = reshape(lineFrame(:,[5,6,11,12]),2,2*2); % 180
            lineL(:,:,4) = reshape(lineFrame(:,[5,6,7,8]),2,2*2); % 270 (L)
            
            lineT(:,:,1) = reshape(lineFrame(:,[1,2,9,10]),2,2*2); % 0 (T)
            lineT(:,:,2) = reshape(lineFrame(:,[3,4,11,12]),2,2*2); % 90
            lineT(:,:,3) = reshape(lineFrame(:,[5,6,9,10]),2,2*2); % 180
            lineT(:,:,4) = reshape(lineFrame(:,[3,4,7,8]),2,2*2); % 270  
              
        end
        %------------------------------
        function drawTargLines(win,targDegree,targEccentricity,lineX,targColor,center,penWidth)
            % targDegree, targEccentricity, penWidth, lineX, and center in pixels
            % targDegree and targEccentricity should have the same size,
            % the length of them indicate the number of lineX
            [targCenter(1,:), targCenter(2,:)] = pol2cart(targDegree,targEccentricity);
            nLineX = size(targCenter,2); % number of LineX as a unit
            nLines = size(lineX,2); % number of single lines in LineX
            lineXs = repmat(lineX,1,nLineX);
            for i = 1:nLineX 
                lineXs(1,1+(i-1)*nLines:i*nLines) = lineX(1,:)+repmat(targCenter(1,i),1,nLines);
                lineXs(2,1+(i-1)*nLines:i*nLines) = lineX(2,:)+repmat(targCenter(2,i),1,nLines);
            end
            Screen('Drawlines', win, lineXs, penWidth, targColor, center); 
        end  
        %-------------------------------------
        function dotXY = cohMotCue(cueDegree,cueDur,dotN,dotSpeed,dotArea,coh)
            % draw coherent motion cue
            % cueDegree in radian
            % cueDur in frames
            % dotSpeed, dotArea in pixels
            % coh in [0,1]
            dotNCoh = round(dotN*coh);
            dotNRan = dotN-round(dotN*coh);
            % initialize dot positions
            theta = -pi+2.*pi.*rand(1,dotN);
            radius = 0+dotArea.*rand(1,dotN);
            dotXY = zeros(2,dotN,cueDur);
            [dotXY(1,:,1),dotXY(2,:,1)] = pol2cart(theta,radius);
            % initialize dot moving direction
            if coh>0 && coh<1
                randMovDir = -pi+2.*pi.*rand(1,dotNRan);
                movDir = [ones(1,dotNCoh).*cueDegree,randMovDir];
            elseif coh==1
                movDir = ones(1,dotN).*cueDegree;
            elseif coh==0
                movDir = -pi+2.*pi.*rand(1,dotN);
            end
            % moving
            for i = 2:cueDur
                % judge whether the dots move outside of the dot area
                % if so, reset the dot position to opposite side
                [movTheta,movR] = cart2pol(dotXY(1,:,i-1),dotXY(2,:,i-1));
                 r_out_edge = find(movR>dotArea);
                 if ~isempty(r_out_edge)
                     movR(r_out_edge) = dotArea;
                     movTheta(r_out_edge) = pi+movTheta(r_out_edge);
                 end
                 [dotXY(1,:,i-1),dotXY(2,:,i-1)] = pol2cart(movTheta,movR);
                 [X,Y] = pol2cart(movDir,dotSpeed); 
                 dotXY(:,:,i) = dotXY(:,:,i-1)+[X;Y];
            end    
        end
        %-------------------------------------
        function [theta,rad] = cueTargCenterTheta(targDiameter,cueTargEcc,choosePos)
            % calcluate possible theta and rads in polar space for cue/target 
            % on the input cueTargetEccentricty and targetDiameter
            % choosePos: the proportion of positions to be chosen (0-1),can
            % be a scalor, or a vector. If it is a vector, the sum should
            % equal 1.
            % two nearby positions at least separated targ Diameter*sqrt(2)
            angleSeparated = targDiameter*sqrt(2)./cueTargEcc; % in radian
            % possible cue and target positions in coordinate space
            nEcc = numel(cueTargEcc);
            nGroup = length(choosePos);
            if nGroup>1 && sum(choosePos)~=1
                fprintf('Warning: the sum of choosePos should equal 1!!!\n');
            end                 
            theta = cell(nGroup,nEcc); % row: subgroup of thetas; column: number of eccentricity layers
            rad = cell(nGroup,nEcc);
            for i = 1:nEcc
                thetaEcc = linspace(0, 2*pi, floor(2.*pi./angleSeparated(i)));
                jitter = 2*pi.*rand; 
                thetaEcc = thetaEcc+jitter; % jitter to increase variability
                posNum = length(thetaEcc);               
                % randomly choose subgroups thetas of cue/targ positions on each ecc
                thetaEcc = Shuffle(thetaEcc);
                for j = 1:nGroup
                    cumsumChoosePos = round(cumsum(choosePos).*posNum);
                    if j == 1
                        theta{j,i} = thetaEcc(1:cumsumChoosePos(j));
                    else
                        theta{j,i} = thetaEcc(cumsumChoosePos(j-1):cumsumChoosePos(j));
                    end
                    rad{j,i} = repmat(cueTargEcc(i),1,length(theta{j,i}));
                end
            end
        end
        %--------------------------------------
        function drawReward(win,scores,barWidth,barHeight,rewardColor)
            % blank tank                
            emptyRect = [0,0,barWidth,barHeight];                           
            Screen('FrameRect', win, rewardColor, emptyRect);
            % filled in reward scores
            rewardRect = [0,0,scores,barHeight];
            Screen('FillRect', win, rewardColor, rewardRect);
        end       
        %-------------------------------------
        function score = rewardCalculate(trialTable,trialID,targDegree,rewardDegree,rewardUnit,rewardRatio,rewardRT)
            % calculate reward score
            RTs = trialTable.RT(1:trialID);
            % remove nan
            RTs = RTs(~isnan(RTs));
            weightRT = 1-normcdf(zscore(RTs),0,1);
            if ~isequal(targDegree,rewardDegree) % referece-reward
                if rewardRT % reward shorter RT
                    score = rewardUnit.*weightRT(end);
                else
                    score = rewardUnit;
                end                                   
            else % reward
                if rewardRT
                    score = rewardUnit.*(1+(rewardRatio-1).*weightRT(end));  
                else
                    score = rewardUnit.*(1+(rewardRatio-1).*rand);
                end
            end        
        end
        %-------------------------------------
        function isRunning = feedbackReward(win,scoreTrial,scoreSum,barWidth,barHeight,center,rewardColor,trialID)
            isRunning = 1; % default
            % feedback reward score to subjects
            if scoreSum<barWidth
                expRun.drawReward(win,scoreSum,barWidth,barHeight,rewardColor); % reward feedback
                performance = sprintf('win: %0.2f',scoreTrial);
                DrawFormattedText(win, performance,'center',center(2)*0.9,rewardColor);
                Screen('Flip', win);   
                WaitSecs(1);
            else
                 textstr='You win! Congratulations!';
                 performance = sprintf('You have spent %d rounds to win the game! Good Job!',trialID);
                 DrawFormattedText(win, textstr,'center','center',200);
                 DrawFormattedText(win, performance, 'center',center(2)-300,200);
                 Screen('Flip', win);
                 WaitSecs(2);
                 isRunning = 0;
                 fprintf('Task finished by Subject! \n');
            end
        end
        %-------------------------------------
        function [resp,rt,isRunning] = keyPressToResp(t0, key, respTimeLimit)
            % initialize
            isRunning = 1; % continue 
            resp = nan;
            rt = nan;                  
            altCheck = 1;
            [keyisdown, ~, keycode] = KbCheck; 
            % wait until response or time ellapse
            while GetSecs-t0<respTimeLimit && (~keyisdown || ...
                 (keyisdown && altCheck && keycode(key.esc)~=1))
                % check key press
                [keyisdown, secs, keycode] = KbCheck;
                for i = 1:length(key.alt) 
                    % check if any key in key.alt{} is pressed
                    altCheck = and(altCheck,keycode(key.alt{i})~=1) ;
                    if keycode(key.alt{i})==1 
                         resp = i;  
                         rt = secs-t0;
                    end
                end
                if keycode(key.esc)==1
                    isRunning = 0;
                    break;
                end
            end
        end
        %-----------------------------------
        function quitPrac = rest(win,trialID,trialTable,blockTrialNum,state,key)
             quitPrac = 0;
             trialNum = height(trialTable);
            if trialID<trialNum && mod(trialID,blockTrialNum)==0
                % break
                Screen('TextSize',win, 30);
                if strcmp(state,'test')
                    textstr='Take a break!';
                    DrawFormattedText(win, textstr,'center','center',200);
                    Screen('Flip', win);
                    % sent feedback to workspace
                    % average ACC and RT
                    avgACC = mean(trialTable.ACC(1:trialID));
                    avgRT  = mean(trialTable.RT(1:trialID));
                    fprintf('Test Feedback: mean ACC:%0.2f, mean RT:%0.2f \n',avgACC,avgRT);
                    WaitSecs(15);
                    textstr='Press Continue key to continue!';
                    DrawFormattedText(win, textstr,'center','center',200);
                    Screen('Flip', win);
                    expRun.keyPressToContinue(key);  
                else
                    textstr='Practice End! \n \n Press Continue key to continue Practice! \n \n Press ESC key to Test \n \n';
                    DrawFormattedText(win, textstr,'center','center',200);
                    % sent feedback to workspace
                    % average ACC and RT
                    avgACC = mean(trialTable.ACC(1:trialID));
                    avgRT  = mean(trialTable.RT(1:trialID));
                    fprintf('Practice Feedback: mean ACC:%0.2f, mean RT:%0.2f \n',avgACC,avgRT);
                    Screen('Flip', win);
                    [~, quitPrac] = expRun.keyPressToContinue(key);                  
                end
            elseif trialID==trialNum 
                % end
                Screen('TextSize',win, 30);
                textstr='The END! \n \n Thank you for participating!';
                DrawFormattedText(win, textstr,'center','center',200);
                Screen('Flip', win);
                WaitSecs(2);
            end 
        end
        %-----------------------------------
        function dataSave(isRunning,expName,subInfo,trialTable,cfg)
            % save data
            if isRunning~=0 && ~strcmp(cfg.state,'practice')
                % make directory if not exist
                savePath = strcat(pwd,strcat('\Data\',expName,'\'));
                if ~isfolder(savePath)
                    mkdir(savePath);
                end
                dataFilename = strcat(subInfo{1},'-',subInfo{3},'-',subInfo{2},'-',datestr(now,30)); 
                save(strcat(savePath,dataFilename),'subInfo','trialTable','cfg');  
                % save table to csv
                writetable(trialTable,strcat(savePath,dataFilename,'.csv'));
            end
        end

    end

end