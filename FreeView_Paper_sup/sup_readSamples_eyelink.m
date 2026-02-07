function M = sup_readSamples_eyelink(smpPath)
% SUP_READSAMPLES_EYELINK  Robust sample reader for EyeLink monocular ASC
% Returns a numeric matrix with 4 columns (monocular only):
%   [time, xp, yp, pupilSize]
% Otherwise raises error.

if ~exist(smpPath,'file')
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

rows = {}; r = 0;
cleanupObj = onCleanup(@() fclose(fid));

% Read sample rows
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    lntrim = strtrim(ln);
    
    % Skip non-sample lines
    if isempty(lntrim)
        continue
    end
    if startsWith(upper(lntrim),'MSG') || startsWith(upper(lntrim),'SAMPLES') || startsWith(upper(lntrim),'EVENTS') || startsWith(upper(lntrim),'START')
        continue
    end
    
    % Sample lines start with number/dot/minus
    if ~(isstrprop(lntrim(1),'digit') || lntrim(1)=='.' || lntrim(1)=='-')
        continue
    end
    
    parts = regexp(lntrim,'\s+','split');
    
    % Expect exactly 4 columns for monocular: TIME X Y PupilSize
    if numel(parts) < 4
        continue
    end
    
    % Parse 4 columns
    nums = nan(1, 4);
    for k = 1:4
        if strcmp(parts{k},'.') || strcmpi(parts{k},'NaN')
            nums(k) = NaN;
        else
            val = str2double(parts{k});
            if ~isnan(val)
                nums(k) = val;
            end
        end
    end
    
    % Require at least TIME (column 1) present
    if isnan(nums(1))
        continue
    end
    
    row = nums(1:4);
    r = r + 1;
    rows{r,1} = row; %#ok<AGROW>
end

if isempty(rows)
    error('No valid sample rows parsed from %s. Expected monocular format: TIME X Y PupilSize', smpPath);
end

M = cell2mat(rows);
end
