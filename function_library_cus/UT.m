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
            % [r, theta; r2, theta2;...] also supported 
            if nargin==2
                distance = obj.default_distance();
            end
            if size(PolarCoor,2)>2 && size(PolarCoor,1)==2
                warning('Be careful! Please check the dim of PolarCoor in UT.Rect2Pol')
                PolarCoor = transpose(PolarCoor);
            end
            x = PolarCoor(:,1).*cosd(PolarCoor(:,2));  % in degree
            y = PolarCoor(:,1).*sind(PolarCoor(:,2));  % in degree
            RectCoor = obj.deg2pix([x,y], distance);
        end

        function PolarCoor = Rect2Pol(obj, RectCoor, distance)
            % [x, y] in pixel to [r, theta] in degree
            % [x1, y1; x2, y2;...] also supported 
            if size(RectCoor,2)>2 && size(RectCoor,1)==2
                warning('Be careful! Please check the dim of RectCoor in UT.Rect2Pol')
                RectCoor = transpose(RectCoor);
            end
            if nargin==2
                distance = obj.default_distance();
            end
            r = sqrt(sum(RectCoor.^2, 2));  % in pixel
            theta = atan2d(RectCoor(:,2),RectCoor(:,1));  % in degree
            PolarCoor = [obj.pix2deg(r, distance),theta];
        end
    end
end

