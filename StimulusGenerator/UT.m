classdef UT % Unit Transformer (单位转换器)
    % 用于视觉刺激实验中的单位转换（像素、厘米、视角度数等）
    % For unit conversion in visual stimulation experiments (pixels, cm, degrees, etc.)
    %
    % 使用示例 / Example usage:
    %   ut = UT(53.5, 1920, 68);  % 屏幕宽53.5cm, 1920像素, 观察距离68cm
    %   pixValue = ut.deg2pix(5);  % 将5度视角转换为像素
    %
    % Yuxin Lu, IPCAS
    % luyx@psych.ac.cn
    
    properties
        width       % 屏幕宽度 (厘米) / Screen width in cm
        Pwidth      % 水平像素数 / Number of horizontal pixels
        distance    % 观察距离 (厘米) / Viewing distance in cm
        ppcm        % 每厘米像素数 / Pixels per cm
        rndPix      % 是否四舍五入像素值 / Whether to round pixel values
    end

    methods
        function obj = UT(width, Pwidth, distance, rndPix)
            % 构造函数 / Constructor
            % width: 屏幕宽度(cm) / Screen width (cm)
            % Pwidth: 水平像素数 / Horizontal pixels
            % distance: 观察距离(cm) / Viewing distance (cm)
            % rndPix: 是否四舍五入(默认true) / Whether to round (default true)
           if nargin >= 3
                obj.distance = distance;
           end
           if nargin == 4
               obj.rndPix = rndPix;
           else
               obj.rndPix = true;
           end
           obj.width = width;
           obj.Pwidth = Pwidth;
           obj.ppcm = Pwidth/width;
        end

        function pix = cm2pix(obj, cm)
            % 厘米转像素 / Convert cm to pixels
            if obj.rndPix
                pix = round(cm * obj.ppcm);
            else
                pix = cm * obj.ppcm;
            end
        end

        function cm = pix2cm(obj, pix)
            % 像素转厘米 / Convert pixels to cm
            cm = pix / obj.ppcm;
        end

        function distance = default_distance(obj)
            if ~isempty(obj.distance)   
                distance = obj.distance;
            else
                error("观察距离尚未设定！/ Distance not set yet!");
            end
        end

        function cm = rad2cm(obj, rad, distance) 
            % 弧度转厘米 / Convert radians to cm
            if nargin == 2
                distance = obj.default_distance();
            end
            cm = tan(rad) * distance;
        end

        function pix = rad2pix(obj, rad, distance)
            % 弧度转像素 / Convert radians to pixels
            if nargin == 2
                distance = obj.default_distance();
            end
            pix = obj.cm2pix(obj.rad2cm(rad, distance));
        end

        function rad = cm2rad(obj, cm, distance)
            % 厘米转弧度 / Convert cm to radians
            if nargin == 2
                distance = obj.default_distance();
            end            
            rad = atan(cm/distance);
        end

        function rad = pix2rad(obj, pix, distance)
            % 像素转弧度 / Convert pixels to radians
            if nargin == 2
                distance = obj.default_distance();
            end
            rad = obj.cm2rad(obj.pix2cm(pix), distance);
        end

        function cm = deg2cm(obj, deg, distance)
            % 度数转厘米 / Convert degrees to cm
            if nargin == 2
                distance = obj.default_distance();
            end
            cm = obj.rad2cm(deg2rad(deg), distance);
        end

        function pix = deg2pix(obj, deg, distance)
            % 度数转像素 / Convert degrees to pixels
            if nargin == 2
                distance = obj.default_distance();
            end
            pix = obj.cm2pix(obj.deg2cm(deg, distance));
        end

        function deg = cm2deg(obj, cm, distance)
            % 厘米转度数 / Convert cm to degrees
            if nargin == 2
                distance = obj.default_distance();
            end
            deg = rad2deg(obj.cm2rad(cm, distance));
        end

        function deg = pix2deg(obj, pix, distance)
            % 像素转度数 / Convert pixels to degrees
            if nargin == 2
                distance = obj.default_distance();
            end
            deg = rad2deg(obj.pix2rad(pix, distance));
        end
        
        function RectCoor = Pol2Rect(obj, PolarCoor, distance)
            % 极坐标转直角坐标 / Polar to Rectangular coordinates
            % 输入: [r, theta] (度数) / Input: [r, theta] (degrees)
            % 输出: [x, y] (像素) / Output: [x, y] (pixels)
            % 支持多行输入 / Supports multiple rows
            if nargin == 2
                distance = obj.default_distance();
            end
            if size(PolarCoor, 2) > 2 && size(PolarCoor, 1) == 2
                warning('请检查PolarCoor的维度！/ Please check PolarCoor dimensions!');
                PolarCoor = transpose(PolarCoor);
            end
            x = PolarCoor(:, 1) .* cosd(PolarCoor(:, 2));
            y = PolarCoor(:, 1) .* sind(PolarCoor(:, 2));
            RectCoor = obj.deg2pix([x, y], distance);
        end

        function PolarCoor = Rect2Pol(obj, RectCoor, distance)
            % 直角坐标转极坐标 / Rectangular to Polar coordinates
            % 输入: [x, y] (像素) / Input: [x, y] (pixels)
            % 输出: [r, theta] (度数) / Output: [r, theta] (degrees)
            % 支持多行输入 / Supports multiple rows
            if size(RectCoor, 2) > 2 && size(RectCoor, 1) == 2
                warning('请检查RectCoor的维度！/ Please check RectCoor dimensions!');
                RectCoor = transpose(RectCoor);
            end
            if nargin == 2
                distance = obj.default_distance();
            end
            r = sqrt(sum(RectCoor.^2, 2));
            theta = atan2d(RectCoor(:, 2), RectCoor(:, 1));
            PolarCoor = [obj.pix2deg(r, distance), theta];
        end
    end
end
