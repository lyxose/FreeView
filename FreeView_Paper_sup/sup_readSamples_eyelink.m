function M = sup_readSamples_eyelink(smpPath)
% SUP_READSAMPLES_EYELINK  Robust sample reader for EyeLink _smp.asc or .asc
% Returns a numeric matrix with at least 7 columns:
%   [time, gxL, gyL, gxR, gyR, pupilL, pupilR]
%
% It tolerates header/comment lines and only keeps lines with >=7 numerics.

if ~exist(smpPath,'file')
    % fallback to base asc
    [p,b,~] = fileparts(smpPath);
    us = strfind(b,'_');
    if ~isempty(us)
        base = b(1:us(1)-1);
    else
        base = b;
    end
    baseAsc = fullfile(p, [base '.asc']);
    if exist(baseAsc,'file')
        smpPath = baseAsc;
    else
        error('Sample ASC file not found: %s', smpPath);
    end
end

fid = fopen(smpPath,'rt');
assert(fid>0, 'Cannot open sample file: %s', smpPath);

% Try to detect field order from a 'FIELDS' header if present
fields = {};
rows = {}; r = 0;
cleanupObj = onCleanup(@() fclose(fid));
while true
    pos = ftell(fid);
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    lntrim = strtrim(ln);
    if startsWith(upper(lntrim),'FIELDS')
        % e.g., 'FIELDS: TIME L_GAZE_X L_GAZE_Y R_GAZE_X R_GAZE_Y L_PUPIL_SIZE R_PUPIL_SIZE'
        toks = regexp(lntrim, 'FIELDS\s*:?(.*)$','tokens','once');
        if ~isempty(toks)
            fields = strsplit(strtrim(toks{1}));
        end
        break
    end
    % Stop scanning headers when we hit first numeric sample
    if ~isempty(lntrim) && (isstrprop(lntrim(1),'digit') || lntrim(1)=='.' || lntrim(1)=='-')
        fseek(fid, pos, -1); % rewind one line
        break
    end
end

% Build a column map (fallback to common default if absent)
function idx = findField(nameList, pat)
    idx = find(~cellfun('isempty', regexp(upper(nameList), pat)), 1);
end

if ~isempty(fields)
    iT  = findField(fields, '^TIME$');
    iLX = findField(fields, '(^|_)L(\w*)GAZE_X$'); if isempty(iLX), iLX = findField(fields,'L_GAZE_X'); end
    iLY = findField(fields, '(^|_)L(\w*)GAZE_Y$'); if isempty(iLY), iLY = findField(fields,'L_GAZE_Y'); end
    iRX = findField(fields, '(^|_)R(\w*)GAZE_X$'); if isempty(iRX), iRX = findField(fields,'R_GAZE_X'); end
    iRY = findField(fields, '(^|_)R(\w*)GAZE_Y$'); if isempty(iRY), iRY = findField(fields,'R_GAZE_Y'); end
    iLP = findField(fields, '(^|_)L(\w*)PUPIL');
    iRP = findField(fields, '(^|_)R(\w*)PUPIL');
else
    % Default EyeLink order (TIME Lx Ly Rx Ry Lp Rp ...)
    iT=1; iLX=2; iLY=3; iRX=4; iRY=5; iLP=6; iRP=7;
end

% Read sample rows
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    lntrim = strtrim(ln);
    % Skip non-sample lines (MSG, event markers, headers)
    if isempty(lntrim)
        continue
    end
    if startsWith(upper(lntrim),'MSG') || startsWith(upper(lntrim),'SAMPLES') || startsWith(upper(lntrim),'EVENTS')
        continue
    end
    % Accept lines beginning with a number, dot, or minus (for negative times in some logs)
    if ~(isstrprop(lntrim(1),'digit') || lntrim(1)=='.' || lntrim(1)=='-')
        continue
    end
    parts = regexp(lntrim,'\s+','split');
    % Map tokens to numeric, '.' -> NaN
    nums = nan(1, max([iT,iLX,iLY,iRX,iRY,iLP,iRP]));
    for k=1:min(numel(parts), numel(nums))
        if strcmp(parts{k},'.') || strcmpi(parts{k},'NaN')
            nums(k) = NaN;
        else
            val = str2double(parts{k});
            if ~isnan(val)
                nums(k) = val;
            end
        end
    end
    % Require at least TIME present
    if isnan(nums(iT))
        continue
    end
    row = [nums(iT), safeIdx(nums,iLX), safeIdx(nums,iLY), safeIdx(nums,iRX), safeIdx(nums,iRY), safeIdx(nums,iLP), safeIdx(nums,iRP)];
    r = r + 1;
    rows{r,1} = row; %#ok<AGROW>
end

if isempty(rows)
    error('No numeric sample rows parsed from %s', smpPath);
end

M = cell2mat(rows);

function v = safeIdx(a,i)
    if i<=numel(a)
        v = a(i);
    else
        v = NaN;
    end
end
end
