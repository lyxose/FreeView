function [times,what,out] = parseMsgs_Form(msgs)

% find starts and ends
iFix   = find(strcmp('FIX ON',msgs(:,2)));
msgs = msgs(iFix(1):end,:);
iFix   = find(strcmp('FIX ON',msgs(:,2)));

starts = regexp(msgs(:,2),'STIM ON: ([\w\\:\.]+)','tokens');
iStart = find(~cellfun(@isempty,starts));
starts = cat(1,starts{iStart}); starts = cat(1,starts{:});

ends   = regexp(msgs(:,2),'STIM OFF: ([\w\\:\.]+)','tokens');
iEnd   = find(~cellfun(@isempty,ends));
ends   = cat(1,ends{iEnd}); ends = cat(1,ends{:});

starts = regexp(msgs(:,2),'STIM ON: ([\w\\:\.]+)','tokens');
iStart = find(~cellfun(@isempty,starts));
starts = cat(1,starts{iStart}); starts = cat(1,starts{:});

ends   = regexp(msgs(:,2),'STIM OFF: ([\w\\:\.]+)','tokens');
iEnd   = find(~cellfun(@isempty,ends));
ends   = cat(1,ends{iEnd}); ends = cat(1,ends{:});

% some check on if file is valid
assert(length(iFix)==length(iEnd))
assert(all(iFix<iEnd))
try 
    assert(all(strcmp(starts,ends)))
catch me
    disp(me)
    disp(starts(1))
    disp(ends(1))
end
    
% get times of starts and ends
times.fix   = cat(1,msgs{iFix});
times.start = cat(1,msgs{iStart});
times.end   = cat(1,msgs{iEnd});

what    = starts;

% collect all messages and their timestamps between each start and end
out = cell(size(iFix));
for p=1:length(iFix)
    out{p} = msgs(iFix(p):iEnd(p),:);
end
