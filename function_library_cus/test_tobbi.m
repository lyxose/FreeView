function test_tobbi(EThndl)
    EThndl.buffer.start('gaze');
    EThndl.sendMessage('Test function',GetSecs)
    WaitSecs(1)
end
