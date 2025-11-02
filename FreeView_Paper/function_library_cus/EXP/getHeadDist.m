function headDist = getHeadDist(gazeData)  
    % return the stream of eye position (row 1: left; row 2: right; row 3: averaged)
    % unit: centimeter 
    headDist = [gazeData.left.gazeOrigin.inUserCoords(3,:); 
                gazeData.right.gazeOrigin.inUserCoords(3,:)]./10; % in cm
    headDist(3,:) = (headDist(1,:) + headDist(2,:)) ./ 2;
end