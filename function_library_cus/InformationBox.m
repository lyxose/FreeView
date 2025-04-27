function [subjID, session, location, subjName, subjGender, subjAge, threshold] = InformationBox
    infoFilePath = './Data/SubjInfo.csv';
    dateNow = str2double(datestr(datetime,'yyyymmdd'));
    exampleInfo = {1, 1, 'IP', 'MingZipinying', 'M', 22, 0, dateNow,''};
    prompt={'Enter Exp. Session:',...
            'Enter location, IP or CAS or Oth:',...
            'Enter Subject Name (QuanPin):',...
            'Enter Subject Gender  [Man: M; Woman: F]:',...
            'Enter Subject Age:',...
            'Enter Central Contrast Threshold: (unknown: 0)',...
            'Enter Date in "yyyymmdd" form:'};
    name='Experimental Information';
    numlines=1;

    % check existance
    if ~exist(infoFilePath, 'file')
        header = {'subjID', 'session', 'location', 'subjName', 'subjGender', 'subjAge', 'threshold', 'THRdate', 'Note'};
        SubjInfo = cell2table(exampleInfo, "VariableNames", header);
        subjID  = input('Enter Subject ID (int):');
        rowIdx = 1;
        defaultanswer = exampleInfo(2:end-1);
    else
        SubjInfo = readtable(infoFilePath);
        SubjInfo.Note = string(SubjInfo.Note);
        subjID   = input(sprintf('Enter Subject ID (next: %.0f):', max(SubjInfo.subjID)+1)); %     
        rowIdx = find(SubjInfo.subjID == subjID,1);
        if ~isempty(rowIdx)
            defaultanswer = table2cell(SubjInfo(rowIdx, 2:end)); 
            defaultanswer{1} = defaultanswer{1}+1; % default as the next session number to prevent data overwrite or comfusing
        else
            rowIdx = height(SubjInfo)+1;
            defaultanswer = exampleInfo(2:end-1);
        end    
    end

    defaultanswer = cellfun(@(x) num2str(x), defaultanswer, 'UniformOutput', false);
    answer   = inputdlg(prompt,name,numlines,defaultanswer);
    session  = str2double(answer{1});  
    location = answer{2};   % IP for Institute of Psychology, CAS for other institutd in CAS, Oth for other subjects
    subjName  = answer{3};
    subjGender= answer{4};
    subjAge   = str2double(answer{5});
    THRdate = str2double(answer{7});
    threshold = str2double(answer{6});  % measured contrast threshold 
%     if ~answer{6}=='0'
%         try
%             threDat = load(sprintf('./Data/Threshold/EYEDat_Sub%.0f_Ses%.0f_%s_%s_%.0f',subjID, session, location, subjName, THRdate));
%         catch me
%             disp(me)
%             disp('NO eye tracker data for this session was found... Be careful when using arbitrary threshold!')
%             files = dir(fullfile('./Data/Threshold'));
%             fileNames = {files.name};
%             matchedFiles = fileNames(contains(fileNames, sprintf('EYEDat_Sub%.0f_', subjID)));
%             threDat = load(sprintf('./Data/Threshold/%s',matchedFiles{end}));
%             if ~contains(matchedFiles{end}, answer{7})
%                 error('NO eye tracker data for this day was found... Do not use the threshold of the previous day!')
%             end
%         end
%     else
%         threDat = [];
%     end
    % update table
    SubjInfo(rowIdx,:) = {subjID, session, location, subjName, subjGender, subjAge, threshold, THRdate, ''};
    writetable(SubjInfo,infoFilePath);
return