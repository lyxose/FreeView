function checkend

[~,~,kc]=KbCheck;
if kc(KbName('escape'))
    ListenChar(0);
    ShowCursor;
    sca;
    disp("To stop the eye tracker: EThndl.buffer.stop('gaze');")
    error('Interprated by user; ”√ªß÷’÷π≥Ã–Ú');
end