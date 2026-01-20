function [subjID, session, location, subjName, subjGender, subjAge] = InformationBox_FV(subjInfo)
% InformationBox_FV  Collect and manage FreeView subject information
% Similar to InformationBox.m but tailored for FreeView experiments
% ALWAYS shows dialog for confirmation, even if info provided
%
% INPUT:
%   subjInfo - struct with optional fields: subjID, session, location,
%              subjName, subjGender, subjAge
%              Can be empty struct or have partial fields
% OUTPUT:
%   subjID, session, location, subjName, subjGender, subjAge

% Parse input struct and extract provided fields (except subjID which is
% always requested separately to follow the original InformationBox logic).
if nargin < 1 || isempty(subjInfo)
    subjInfo = struct();
end

session_in   = getfield_safe(subjInfo, 'session', []); %#ok<GFLD>
location_in  = getfield_safe(subjInfo, 'location', '');
subjName_in  = getfield_safe(subjInfo, 'subjName', '');
subjGender_in= getfield_safe(subjInfo, 'subjGender', '');
subjAge_in   = getfield_safe(subjInfo, 'subjAge', NaN);

dateNow = str2double(datestr(datetime,'yyyymmdd'));
exampleInfo = {1, 1, 'IP', 'MingZipinying', 'M', 22, dateNow, ''};
infoFilePath = './Data/SubjInfo_FV.csv';
dateNow = str2double(datestr(datetime,'yyyymmdd'));
exampleInfo = {1, 1, 'IP', 'MingZipinying', 'M', 22, dateNow, ''};

% Dialog labels
promptID = {'Enter Subject ID (FV only):'};
promptAll = {
    'Enter Exp. Session:' ...
    'Enter location (IP/CAS/Oth):' ...
    'Enter Subject Name (QuanPin):' ...
    'Enter Subject Gender [M/F]:' ...
    'Enter Subject Age:' ...
    'Enter Date in "yyyymmdd" form:'
};
name = 'FreeView Subject Information';
numlines = 1;

% Initialize or load existing SubjInfo_FV.csv, then collect subjID FIRST
if ~exist(infoFilePath, 'file')
    header = {'subjID', 'session', 'location', 'subjName', 'subjGender', 'subjAge', 'date', 'Note'};
    SubjInfo = cell2table(cell(0, numel(header)), 'VariableNames', header);
    suggestedID = 1;
    % Defaults when no history exists
    defSession  = fallback_numeric(session_in, 1);
    defLocation = fallback_text(location_in, 'IP');
    defName     = fallback_text(subjName_in, '');
    defGender   = fallback_text(subjGender_in, '');
    defAge      = fallback_numeric(subjAge_in, 22);
else
    SubjInfo = readtable(infoFilePath);
    SubjInfo.Note = string(SubjInfo.Note);
    suggestedID = max(SubjInfo.subjID) + 1;
end

% Always ask for FV subjID (separate dialog to mirror original InformationBox)
idDefault = {num2str(suggestedID)};
input('Press Enter to Continue:');
idAns = inputdlg(promptID, name, numlines, idDefault);
if isempty(idAns)
    error('User cancelled subject information input');
end
subjID = str2double(idAns{1});
if isnan(subjID)
    error('Invalid Subject ID');
end

% Build defaults for the full dialog based on table contents
rowIdx = find(SubjInfo.subjID == subjID, 1);
if ~isempty(rowIdx)
    existingRow = table2cell(SubjInfo(rowIdx, 1:end-1));
    if isempty(session_in)
        suggestedSession = existingRow{2} + 1;
    else
        suggestedSession = session_in;
    end
    defaultanswer = {
        num2str(suggestedSession) ...
        char(existingRow{3}) ...
        char(existingRow{4}) ...
        char(existingRow{5}) ...
        num2str(existingRow{6}) ...
        num2str(existingRow{7})
    };
else
    rowIdx = height(SubjInfo) + 1;
    defSession  = fallback_numeric(session_in, 1);
    defLocation = fallback_text(location_in, 'IP');
    defName     = fallback_text(subjName_in, '');
    defGender   = fallback_text(subjGender_in, '');
    defAge      = fallback_numeric(subjAge_in, 22);
    defaultanswer = {
        num2str(defSession) ...
        defLocation ...
        defName ...
        defGender ...
        num2str(defAge) ...
        num2str(dateNow)
    };
end

% NOW show input dialog for remaining fields
answer = inputdlg(promptAll, name, numlines, defaultanswer);
if isempty(answer)
    error('User cancelled subject information input');
end

% Parse answers
session   = str2double(answer{1});
location  = answer{2};
subjName  = answer{3};
subjGender= answer{4};
subjAge   = str2double(answer{5});
dateInput = str2double(answer{6});

% Update/create table row
SubjInfo(rowIdx, :) = {subjID, session, location, subjName, subjGender, subjAge, dateInput, ''};

% Ensure Data directory exists
if ~exist('./Data', 'dir')
    mkdir('./Data');
end

% Write to CSV
writetable(SubjInfo, infoFilePath);

end

% ---- Helpers ----
function val = getfield_safe(s, fname, defaultVal)
if isstruct(s) && isfield(s, fname) && ~isempty(s.(fname))
    val = s.(fname);
else
    val = defaultVal;
end
end

function val = fallback_numeric(v, defaultVal)
if isempty(v) || isnan(v)
    val = defaultVal;
else
    val = v;
end
end

function val = fallback_text(v, defaultVal)
if isempty(v)
    val = defaultVal;
else
    val = v;
end
end
