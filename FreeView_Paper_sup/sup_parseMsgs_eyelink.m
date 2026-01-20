function [times, what, out] = sup_parseMsgs_eyelink(msgs)
% SUP_PARSEMSGS_EYELINK  Parse EyeLink MSG array into trial-aligned times
%   msgs: Nx2 cell {timestamp(ms), text}
% Returns:
%   times.fix   : fix-on timestamps (vector)
%   times.start : stimulus-on timestamps (vector)
%   times.end   : stimulus-off timestamps (vector)
%   what        : cellstr payload after 'STIM ON:' per trial (may be empty)
%   out         : cell array of sub-messages from fix..end per trial

% Extract indices
msgText = msgs(:,2);
msgTime = cell2mat(msgs(:,1));

isFix   = startsWith(msgText,'FIX ON');
isStimOn  = startsWith(msgText,'STIM ON');
isStimOff = startsWith(msgText,'STIM OFF');

idxFix   = find(isFix);
idxOn    = find(isStimOn);
idxOff   = find(isStimOff);

% Walk through messages in order to pair FIX -> STIM ON -> STIM OFF
tri_fix = [];
tri_on  = [];
tri_off = [];
tri_what= {};
tri_out = {};

kFix = 1; kOn = 1; kOff = 1; n = numel(msgText);

while kFix <= numel(idxFix)
    iFix = idxFix(kFix);
    % find first STIM ON after iFix
    iOnCand = idxOn(kOn:end);
    iOnCand = iOnCand(iOnCand > iFix);
    if isempty(iOnCand)
        break
    end
    iOn = iOnCand(1);
    % find first STIM OFF after iOn
    iOffCand = idxOff(kOff:end);
    iOffCand = iOffCand(iOffCand > iOn);
    if isempty(iOffCand)
        break
    end
    iOff = iOffCand(1);

    tri_fix(end+1,1) = msgTime(iFix); %#ok<AGROW>
    tri_on (end+1,1) = msgTime(iOn);  %#ok<AGROW>
    tri_off(end+1,1) = msgTime(iOff); %#ok<AGROW>

    % Extract payload after 'STIM ON:' for this trial (optional)
    m = msgText{iOn};
    tok = regexp(m,'STIM ON:\s*(.*)$','tokens','once');
    if isempty(tok), tok = {''}; end
    tri_what{end+1,1} = tok{1}; %#ok<AGROW>

    % Collect sub-messages from fix..off
    tri_out{end+1,1} = msgs(iFix:iOff,:); %#ok<AGROW>

    % Advance pointers past these indices
    kFix = find(idxFix > iOff, 1, 'first'); if isempty(kFix), kFix = numel(idxFix)+1; end
    kOn  = find(idxOn  > iOff, 1, 'first'); if isempty(kOn ), kOn  = numel(idxOn) +1; end
    kOff = find(idxOff > iOff, 1, 'first'); if isempty(kOff), kOff = numel(idxOff)+1; end
end

% Pack outputs
times.fix   = tri_fix;
times.start = tri_on;
times.end   = tri_off;
what        = tri_what;
out         = tri_out;
end
