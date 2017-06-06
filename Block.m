% Block class
% Software implementation copyright Adrian Harwood 2017

classdef Block < handle
    % Class to hold information about a block. May be extended to represent 
    % patterns for specific items of clothing.

    properties (Access = public)
        % Name of the block
        name = 'default'
    end
    
    properties (Access = private)
        % Privately stored keypoints
        keypointsX = [];
        keypointsY = [];
    end
    
    %% Public methods
    methods (Access = public)
        
        %% Constructor
        function obj = Block(origObj)
            if nargin == 0
                % Default constructor
            elseif nargin == 1
                % Copy constructor
                props = properties(origObj);
                for i = 1 : length(props)
                    obj.(props{i}) = origObj.(props{i});
                end
            end
        end
        
        %% addKeypoint
        function addKeypoint(obj, xy)
            % Add a point to the end of the keypoints list
            obj.keypointsX = [obj.keypointsX; xy(1)];  
            obj.keypointsY = [obj.keypointsY; xy(2)];
        end
        
        %% addKeypointNextTo
        function addKeypointNextTo(obj, xy, adjacent, place)
            % Add a point before or after the adjacent point
                
            % Loop through and find point given
            for i = 1 : length(obj.keypointsX)
                if (obj.keypointsX(i) == adjacent(1) && ...
                        obj.keypointsY(i) == adjacent(2))
                    break;
                end
                % If reached the end then it will be added to the end
            end                

            if (strcmp(place,'before')) 
                i = i;
            elseif (strcmp(place,'after'))
                i = i + 1;
            else
                error('bad argument calling addKeypoint');                    
            end
            
            obj.keypointsX(1:i-1) = obj.keypointsX(1:i-1);
            tmp = obj.keypointsX(i:end);
            obj.keypointsX(i) = xy(1);
            obj.keypointsX(i+1:end+1) = tmp;

            obj.keypointsY(1:i-1) = obj.keypointsY(1:i-1);
            tmp = obj.keypointsY(i:end);
            obj.keypointsY(i) = xy(2);
            obj.keypointsY(i+1:end+1) = tmp;                
        end
        
        %% addDart
        function dartEdges = addDart(obj, line_start, line_end, ...
                dimensionless_position, width, length)
            % Add a dart given the end points of the line one which to add 
            % it, the location along the line, the dart width and length.
            
            % Find equation of line
            direction = line_end - line_start;
            normal = [-direction(2) direction(1)];
            
            % Normalise the direction vectors
            normal  = normal / norm(normal);
            direction = direction / norm(direction);
            
            % Find dart point
            side_length = norm(line_end - line_start);
            point = line_start + (dimensionless_position * side_length)...
                * direction;
            
            % Add the keypoints associated with the dart
            pt1 = point + (width / 2) * direction;
            pt2 = point + length * normal;
            pt3 = point - (width / 2) * direction;
            dartEdges = [pt1; pt2; pt3];
            obj.addKeypointNextTo(pt1, line_end, 'before');
            obj.addKeypointNextTo(pt2, pt1, 'before');
            obj.addKeypointNextTo(pt3, pt2, 'before');
        end        
        
        %% addRightAngleCurve
        function addRightAngleCurve(obj, point1, point2, dir1, dir2)
            % Add a curve given the direction of the two lines which should
            % be connected to at a right angle.
            
            % Flag to indicate we switch coordinates to handle high
            % gradients
            bSwitched = false;
        
            % Normalise direction vectors
            len1 = norm(dir1);
            len2 = norm(dir2);
            dir1 = dir1 / len1;
            dir2 = dir2 / len2;
            
            % Find normals at each point (gradients)
            norm1 = [-dir1(2) dir1(1)];
            dydx1 = -norm1(2) / norm1(1);
            norm2 = [-dir2(2) dir2(1)];
            dydx2 = norm2(2) / norm2(1);
            
            % If gradients are too large then rotate coordinates
            if (abs(dydx1) > 1e6 || abs(dydx2) > 1e6)
                point1 = fliplr(point1);
                point2 = fliplr(point2);
                dir1 = fliplr(dir1);
                dir2 = fliplr(dir2);
                norm1 = [-dir1(2) dir1(1)];
                dydx1 = norm1(2) / norm1(1);
                norm2 = [-dir2(2) dir2(1)];
                dydx2 = norm2(2) / norm2(1);
                bSwitched = true;
            end
            
            % Find coefficients for the cubic spline
            % ax^3 + bx^2 + cx + d
            A = [...
                point1(1)^3 point1(1)^2 point1(1) 1; ...
                point2(1)^3 point2(1)^2 point2(1) 1; ...
                3*point1(1)^2 2*point1(1) 1 0; ...
                3*point2(1)^2 2*point2(1) 1 0];
            b = [point1(2); point2(2); dydx1; dydx2];
            coeff = A \ b;
            
            % Discretise to roughly 2 points per cm
            numPts = ceil(norm(point2 - point1) * 10);
            
            if (bSwitched)
                xpts = linspace(point1(1), point2(1), numPts);
                ypts = coeff(1) * xpts.^3 + coeff(2) * xpts.^2 + coeff(3) * xpts + coeff(4);
                tmp = xpts;
                xpts = ypts;
                ypts = tmp;
                point1 = fliplr(point1);
            else
                xpts = linspace(point1(1), point2(1), numPts);
                ypts = coeff(1) * xpts.^3 + coeff(2) * xpts.^2 + coeff(3) * xpts + coeff(4);
            end
            %figure, plot(xpts, ypts, 'x-'), axis equal;
            
            % Add keypoints           
            xpts = flipud(xpts);
            ypts = flipud(ypts);
            obj.addKeypointNextTo([xpts(2) ypts(2)], ...
                    [point1(1) point1(2)],'after');
            for i = 3 : length(xpts)-1
                obj.addKeypointNextTo([xpts(i) ypts(i)], ...
                    [xpts(i-1) ypts(i-1)],'after');
            end   
            
        end
        
        %% addCircularCurve
        function addCircularCurve(obj, xy1, xy2, h)
            % Add a curve given the height of curve above the centre of a 
            % line joining the two points. Assumes the curve is cut from a 
            % circle.
            
            % Get equation of the line between and hence midpoint
            direction = (xy2 - xy1);
            norm_dir = [-direction(2) direction(1)];
            norm_dir = norm_dir / norm(norm_dir);
            
            % Use normalised direction since the 0.5 is normalised
            midpoint = xy1 + 0.5 * direction;            
            
            % Find third point (use unit vector normal since h in cm)
            xy3 = midpoint + h * norm_dir;
            
            % Solve for missing coefficients
            lam1 = -xy1(1)^2 - xy1(2)^2;
            lam2 = -xy2(1)^2 - xy2(2)^2;
            lam3 = -xy3(1)^2 - xy3(2)^2;
            lam23 = lam2 - lam3;
            lam13 = lam1 - lam3;
            y23 = xy2(2) - xy3(2);
            y13 = xy1(2) - xy3(2);
            x23 = xy2(1) - xy3(1);
            x13 = xy1(1) - xy3(1);
            alp = (lam23 - (y23/y13) * lam13) / ((x23 * y13 - y23 * x13) / y13);
            bet = (lam13 - x13 * alp) / y13;
            c = lam1 - xy1(2) * bet - xy1(1) * alp;
            
            % Convert to standard form of equation for circle
            a = -alp / 2;
            b = -bet / 2;
            radius = sqrt(a^2 + b^2 - c);         
            
            % Discretise equation of circle using 2 points per cm
            th1 = acos((xy1(1) - a) / radius);
            th2 = acos((xy2(1) - a) / radius);
            dcircum = abs(th2 - th1) * radius;
            numPts = ceil(dcircum * 2);
            
            % Specify in polar coordinates
            th = linspace(th1, th2, numPts);
            r = radius * ones(1, length(th));
            [xpts, ypts] = pol2cart(th, r);
            xpts = xpts + a;
            ypts = ypts + b;
            %plot(xpts,ypts,'x-'), axis equal;
            
            % Add keypoints
            xpts = flipud(xpts);
            ypts = flipud(ypts);
            obj.addKeypointNextTo([xpts(2) ypts(2)], ...
                    [xy1(1) xy1(2)],'before');
            for i = 3 : length(xpts)-1
                obj.addKeypointNextTo([xpts(i) ypts(i)], ...
                    [xpts(i-1) ypts(i-1)],'before');
            end
        
        end
        
        %% plotBlock
        function plotBlock(obj)
            % Plots the block to a new figure
            
            figure;
            % Add last first point to last to close contour
            plot([obj.keypointsX; obj.keypointsX(1)],...
                [obj.keypointsY; obj.keypointsY(1)],'-x');
            axis equal;
        end
    end
    
end    