
function stimulus = genStim(winRect, ut, bgContrast, tgContrast, tgCenter, GaborSF, GaborWidth, GaborOrient, bgWidth, seed)
    % to generate gray matrix of stimulus (Gabor in a round 1/f pinkNoise texture)
    % tgCenter (i.e. target center coor.) is in pixel unit.
    % all last 4 parameters are in degree! (GaborSF, GaborWidth, GaborOrient, bgWidth)

    scWidth  = winRect(3);  
    scHeight = winRect(4);
    bgCenter = [scWidth/2, scHeight/2];
    
    background = tPinkNoise(scWidth, seed, bgContrast); % full screen background texture
    background = background(1:scHeight,1:scWidth);
    lambda = ut.deg2pix(1/GaborSF);
    if lambda ~= 0
        Texture = grating(size(background), bgCenter+[tgCenter(1),-tgCenter(2)], ...
                            1/lambda, GaborOrient, tgContrast);
        Texture = winOverlap(background, Texture, ut.deg2pix(GaborWidth), ...
                              bgCenter+[tgCenter(1),-tgCenter(2)], 'cos'); 
    else
        fprintf('Spatial frequency was too large (%.4f) that grating cannot be generated at current screen resolution',GaborSF)
        Texture = background;
    end
    stimulus = winOverlap(zeros([scHeight,scWidth])+0.5, Texture, ...
                          ut.deg2pix(bgWidth), bgCenter, 'hard');
end     