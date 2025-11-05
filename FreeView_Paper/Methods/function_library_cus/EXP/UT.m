classdef UT % Unit Transformer
    properties
        width       % width of screen, in centimeter
        Pwidth      % number of horizontal pixels
        distance    % distance between face and screen, in centimeter, empty if dynamic
        ppcm        % pixels per centimeter
        rndPix      % if true, return the pixel in int
    end

    methods
        function obj = UT(width, Pwidth, distance, rndPix)
           if nargin>=3
                obj.distance = distance;
           end
           if nargin ==4
               obj.rndPix = rndPix;
           else
               obj.rndPix = true;
           end
           obj.width = width;       % width of screen, in centimeter
           obj.Pwidth = Pwidth;     % number of horizontal pixels
           obj.ppcm = Pwidth/width; % pixels per centimeter, approximated by central ppcm
        end

        function pix = cm2pix(obj, cm)
            if obj.rndPix
                pix = round(cm * obj.ppcm); % not invertible at low resolution 
            else
                pix = cm * obj.ppcm;
            end
        end

        function cm = pix2cm(obj, pix)
            cm = pix / obj.ppcm;
        end

        function distance = default_distance(obj)
            if ~isempty(obj.distance)   
                distance=obj.distance;
            else
                error("Distance between screen and face have not been given yet!")
            end
        end

        function cm = rad2cm(obj, rad, distance) 
            if nargin==2
                distance = obj.default_distance();
            end
            cm = tan(rad)*distance;
        end

        function pix = rad2pix(obj, rad, distance)
            if nargin==2
                distance = obj.default_distance();
            end
            pix = obj.cm2pix(obj.rad2cm(rad,distance));
        end

        function rad = cm2rad(obj, cm, distance)
            if nargin==2
                distance = obj.default_distance();
            end            
            rad = atan(cm/distance);
        end

        function rad = pix2rad(obj, pix, distance)
            if nargin==2
                distance = obj.default_distance();
            end
            rad = obj.cm2rad(obj.pix2cm(pix), distance);
        end

        function cm = deg2cm(obj, deg, distance)
            if nargin==2
                distance = obj.default_distance();
            end
            cm = obj.rad2cm(deg2rad(deg), distance);
        end

        function pix = deg2pix(obj, deg, distance)
            if nargin==2
                distance = obj.default_distance();
            end
            pix = obj.cm2pix(obj.deg2cm(deg, distance));
        end

        function deg = cm2deg(obj, cm, distance)
            if nargin==2
                distance = obj.default_distance();
            end
            deg = rad2deg(obj.cm2rad(cm), distance);
        end

        function deg = pix2deg(obj, pix, distance)
            if nargin==2
                distance = obj.default_distance();
            end
            deg = rad2deg(obj.pix2rad(pix, distance));
        end
        
        function RectCoor = Pol2Rect(obj, PolarCoor, distance)
            % [r, theta] in degree to [x, y] in pixel
            if nargin==2
                distance = obj.default_distance();
            end
            x = PolarCoor(1)*cosd(PolarCoor(2));  % in degree
            y = PolarCoor(1)*sind(PolarCoor(2));  % in degree
            RectCoor = obj.deg2pix([x,y], distance);
        end

        function PolarCoor = Rect2Pol(obj, RectCoor, distance)
            % [x, y] in pixel to [r, theta] in degree
            if nargin==2
                distance = obj.default_distance();
            end
            r = norm(RectCoor);  % in pixel
            if RectCoor(1)==0
                theta = 180 - 90*sign(RectCoor(2)); % 90 or 270
            else
                theta = atand(RectCoor(2)/RectCoor(1));  % in degree
            end
            if RectCoor(1)<0
                theta = theta+180;
            end
            PolarCoor = [obj.pix2deg(r, distance),theta];
        end
    end
end

