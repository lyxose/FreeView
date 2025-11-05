function drawcenteredtext_dot(w,str,x,y,color,textsize)
% Author: Zijian Chen @ SJTU
str=double(str);
oldsize=Screen('TextSize',w,textsize);
drect=Screen('TextBounds',w,str);

sx=x-drect(3)/2;
sy=y-drect(4)/2;
Screen('DrawText',w,str,sx,sy,color);
Screen('TextSize',w,oldsize);