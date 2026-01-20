function msgs = sup_loadMsgs_eyelink(msgPath)
% SUP_LOADMSGS_EYELINK  Load EyeLink ASC messages into Nx2 cell array
%   msgs(:,1) = timestamp (double, ms)
%   msgs(:,2) = message string
%
% It accepts either a dedicated *_msg.asc file or a single *.asc file.
% Lines are expected in edf2asc format like:
%   MSG 1234567 STIM ON: ...
%   MSG 1237890 FIX ON
%
% If msgPath does not exist, tries a fallback base ".asc" path.

if ~exist(msgPath,'file')
    [p,b,~] = fileparts(msgPath);
    baseAsc = fullfile(p, [b(1:max(1,strfind(b,'_')-1)) '.asc']);
    if exist(baseAsc,'file')
        msgPath = baseAsc;
    else
        error('Message ASC file not found: %s', msgPath);
    end
end

fid = fopen(msgPath,'rt');
assert(fid>0, 'Cannot open message file: %s', msgPath);
C = {};
try
    % Stream line-by-line to be robust to mixed content
    i = 0;
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        lntrim = strtrim(ln);
        % EyeLink messages start with 'MSG'
        if startsWith(upper(lntrim),'MSG')
            % Format: MSG <time> <message...>
            parts = regexp(lntrim,'^MSG\s+(\d+)\s+(.*)$','tokens','once');
            if ~isempty(parts)
                i = i+1;
                t = str2double(parts{1});
                m = strtrim(parts{2});
                C{i,1} = t; %#ok<AGROW>
                C{i,2} = m; %#ok<AGROW>
            end
        end
    end
finally
    fclose(fid);
end

if isempty(C)
    warning('No MSG lines found in %s', msgPath);
end
msgs = C;
end
