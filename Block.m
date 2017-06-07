% Block class
% Class represents a block defined as a series of connected keypoints
% representing the geometrical design.
% Software implementation copyright Adrian Harwood 2017.
% The University of Manchester, UK.

classdef Block < handle
% Class to hold information about a block. May be extended to represent 
% patterns for specific items of clothing.

    properties (Access = public)
        % Name of the block.
        name = 'default'
    end
    
    properties (Access = private)
        % Privately stored keypoints ordered in anti-clockwise sense for
        % correct plotting.
        keypointsX = [];
        keypointsY = [];
    end
    
    %% Public methods
    methods (Access = public)
        
        %% Constructor
        function this = Block(origObj)
            if nargin == 0
                % Default constructor.
                
            elseif nargin == 1
                % Copy constructor.
                
                % The line below only copies public data which is not
                % useful to us here as we want to clone the whole object.
                %props = properties(origObj);
                
                % This works as long as this class has no handles to other
                % handle class in which case the handles would still point
                % to the old object. Also MATLAB complains that we are
                % trying to copy private data. Should be a way to suppress
                % this warning.
                props = fieldnames(struct(origObj));
                
                % Copy the properties one at a time into the new copy
                for i = 1 : length(props)
                    this.(props{i}) = origObj.(props{i});
                end
            end
        end
        
        %% getKeypointNumber
        function ptnum = getKeypointNumber(this, xy)
        % Method to find which keypoint number corresponds to the position
        % given (within floating point tolerances)
        
            % Loop through and find point given (to eps limit).
            i = 1;
            while (i <= length(this.keypointsX))
                if (...
                        abs(this.keypointsX(i) - xy(1)) < eps && ...
                        abs(this.keypointsY(i) - xy(2)) < eps ...
                    )
                    % Point found.
                    break;
                end
                
                % Increment counter
                i = i + 1;
            end 
            
            % If i = length + 1, point not found so return -1
            if (i > length(this.keypointsX))
                ptnum = -1;
            else
                % Point found and position = i
                ptnum = i;
            end            
        end
        
        %% addKeypoint
        function addKeypoint(this, xy)
        % Append a point to the end of the keypoints list.
            this.keypointsX = [this.keypointsX; xy(1)];  
            this.keypointsY = [this.keypointsY; xy(2)];
        end
        
        %% addKeypointNextTo
        function addKeypointNextTo(this, xy, adjacent, place)
        % Add a point before or after the adjacent point.
                
            % Get point number of adjacent point.
            i = this.getKeypointNumber(adjacent);
            
            % If point not found then error.
            if (i == -1)
                error('Adjacent point not found in keypoints list');
            end

            % Correct position in list if adding after.
            if (strcmp(place,'after'))
                i = i + 1;
            elseif (~strcmp(place,'before'))
                error('Bad string argument calling addKeypointNextTo');
            end
            
            % Split points list and insert point.
            this.keypointsX(1:i-1) = this.keypointsX(1:i-1);
            tmp = this.keypointsX(i:end);
            this.keypointsX(i) = xy(1);
            this.keypointsX(i+1:end+1) = tmp;

            this.keypointsY(1:i-1) = this.keypointsY(1:i-1);
            tmp = this.keypointsY(i:end);
            this.keypointsY(i) = xy(2);
            this.keypointsY(i+1:end+1) = tmp;                
        end
        
        %% addDart
        function dartEdges = addDart(this, line_start, line_end, ...
                dimensionless_position, width, length)
        % Add a dart given the end points of the line on which to add 
        % it, the location along the line, the dart width and length.
        % Start and end points must be specified in the strict 
        % anti-clockwise order for connectivity to be correct for plotting.
            
            % Find equation of line.
            direction = line_end - line_start;
            normal = [-direction(2) direction(1)];
            
            % Normalise the direction vectors.
            normal  = normal / norm(normal);
            direction = direction / norm(direction);
            
            % Find dart point.
            side_length = norm(line_end - line_start);
            point = line_start + (dimensionless_position * side_length)...
                * direction;
            
            % Add the keypoints associated with the dart.
            pt1 = point - (width / 2) * direction;
            pt2 = point + length * normal;
            pt3 = point + (width / 2) * direction;
            dartEdges = [pt1; pt2; pt3];
            this.addKeypointNextTo(pt1, line_start, 'after');
            this.addKeypointNextTo(pt2, pt1, 'after');
            this.addKeypointNextTo(pt3, pt2, 'after');
        end        
        
        %% addRightAngleCurve
        function addRightAngleCurve(this, startpt, endpt, dirstart, dirend, dirnorm)
        % Add a curve between the start and end points which meets each
        % bounding line at 90 degrees. Direction of both bounding lines
        % must be given as well as whether normal should be to the right or
        % left for each bounding edge. This is indiciated as a true or 
        % false 'dirnorm' flag. Note: the normals need to have a magintude
        % that makes sense in the context of the start and end point
        % specification order as well as the world space as this is the
        % system used to specify the equation of the curve. More often than
        % not, the normals will be pointing in the same direction so the
        % flags are usually set to the same value.
        
        
            % To ensure the process is more robust we first map the two end
            % points and their side directions onto a reference X axis.
            % This is done by first shifting the start point to the
            % reference origin then rotating the points and directions such
            % that the start bounding line is coincident with the reference
            % Y axis.

            % Shift the start and end points
            refstart = startpt - startpt;
            refend = endpt - startpt;

            % Find rotation angle
            rotang = acos(...
                dot(dirstart, [0 1]) / (norm(dirstart) * norm([0 1])) ...
                );

            % Construct rotation matrix
            R = [cos(rotang) -sin(rotang); sin(rotang) cos(rotang)];

            % Rotate the direction vectors and position vectors
            refdirstart = R * dirstart';
            refdirend = R * dirend';
            refstart = R * refstart';
            refend = R * refend';            

            % Compute the normal vectors in reference system
            if (dirnorm(1))
                % Right hand normal.
                refnormstart = [refdirstart(2) -refdirstart(1)];
            else
                % Left hand normal.
                refnormstart = [-refdirstart(2) refdirstart(1)];
            end

            if (dirnorm(2))
                refnormend = [refdirend(2) -refdirend(1)];
            else
                refnormend = [-refdirend(2) refdirend(1)];
            end

            % Compute the required gradients of the curve.
            refdydxstart = refnormstart(2) / refnormstart(1);
            refdydxend = refnormend(2) / refnormend(1);

            % Find coefficients for the cubic spline:
            % ax^3 + bx^2 + cx + d
            % By using two end point conditions and two gradient conditions
            % to define a set of simultaneous equations.
            A = [...
                refstart(1)^3 refstart(1)^2 refstart(1) 1; ...
                refend(1)^3 refend(1)^2 refend(1) 1; ...
                3*refstart(1)^2 2*refstart(1) 1 0; ...
                3*refend(1)^2 2*refend(1) 1 0];
            b = [refstart(2); refend(2); refdydxstart; refdydxend];

            % Solve system to get coefficients of spline.
            coeff = A \ b;

            % Discretise to roughly 1 point per cm assuming a straight
            % line between points for simplicity.
            numPts = ceil(norm(refend - refstart) * 1);

            % Find points on curve by seeding x.
            % This might not always be robust -- should use a local
            % curvilinear coordinate system really.
            xpts = linspace(refstart(1), refend(1), numPts);
            ypts = coeff(1) * xpts.^3 + coeff(2) * xpts.^2 + coeff(3) * xpts + coeff(4);
            %figure, subplot(3,1,1), plot(xpts, ypts, 'x-'), axis equal;
            
            % Reverse rotation
            Ri = [cos(-rotang) -sin(-rotang); sin(-rotang) cos(-rotang)];
            for i = 1 : numPts
                tmp = [xpts(i); ypts(i)];
                tmp = Ri * tmp;
                xpts(i) = tmp(1);
                ypts(i) = tmp(2);
            end
            %subplot(3,1,2), plot(xpts, ypts, 'x-'), axis equal;
            
            % Reverse shift
            xpts = xpts + startpt(1);
            ypts = ypts + startpt(2);
            %subplot(3,1,3), plot(xpts, ypts, 'x-'), axis equal;
            
            % Add keypoints           
            this.addKeypointNextTo([xpts(2) ypts(2)], ...
                    [startpt(1) startpt(2)],'after');
            for i = 3 : length(xpts)-1
                this.addKeypointNextTo([xpts(i) ypts(i)], ...
                    [xpts(i-1) ypts(i-1)],'after');
            end   
            
        end
        
        %% addCircularCurve
        function addCircularCurve(this, startpt, endpt, h, dirnorm)
        % Add a curve given the height of curve above the centre of a 
        % line joining the two points. Assumes the curve is cut from a 
        % circle and hence given points are on the circle circumference.
        % Direction of normal is indicated by boolean value -- true for
        % right hand normal and false for left hand normal -- which tells
        % the method which way to curve. Start and end points must be
        % specified in the strict anti-clockwise order.
            
            % Get equation of the line between and hence midpoint
            direction = (endpt - startpt);
            if (dirnorm)
                % Right hand normal
                norm_dir = [direction(2) -direction(1)];
            else
                % Left hand normal
                norm_dir = [-direction(2) direction(1)];
            end
            norm_dir = norm_dir / norm(norm_dir);
            
            % Use normalised direction since the 0.5 is normalised
            midpt = startpt + 0.5 * direction;            
            
            % Find third point (use unit vector normal since h in cm)
            crestpt = midpt + h * norm_dir;
            
            % Solve for missing coefficients
            lam1 = -startpt(1)^2 - startpt(2)^2;
            lam2 = -endpt(1)^2 - endpt(2)^2;
            lam3 = -crestpt(1)^2 - crestpt(2)^2;
            lam23 = lam2 - lam3;
            lam13 = lam1 - lam3;
            y23 = endpt(2) - crestpt(2);
            y13 = startpt(2) - crestpt(2);
            x23 = endpt(1) - crestpt(1);
            x13 = startpt(1) - crestpt(1);
            alp = (lam23 - (y23/y13) * lam13) / ((x23 * y13 - y23 * x13) / y13);
            bet = (lam13 - x13 * alp) / y13;
            gam = lam1 - startpt(2) * bet - startpt(1) * alp;
            
            % Convert to standard form of equation for circle
            centrex = -alp / 2;
            centrey = -bet / 2;
            radius = sqrt(centrex^2 + centrey^2 - gam);         
            
            % Discretise equation of circle using 1 points per cm
            th1 = acos((startpt(1) - centrex) / radius);
            th2 = acos((endpt(1) - centrex) / radius);
            dcircum = abs(th2 - th1) * radius;
            numPts = ceil(dcircum * 1);
            
            % Specify in polar coordinates
            th = linspace(th1, th2, numPts);
            r = radius * ones(1, length(th));
            [xpts, ypts] = pol2cart(th, r);
            xpts = xpts + centrex;
            ypts = ypts + centrey;
            %plot(xpts,ypts,'x-'), axis equal;
            
            % Add keypoints (skip first point, assume we already added it)
            this.addKeypointNextTo([xpts(2) ypts(2)], ...
                    [startpt(1) startpt(2)],'after');
            for i = 3 : length(xpts)-1
                this.addKeypointNextTo([xpts(i) ypts(i)], ...
                    [xpts(i-1) ypts(i-1)],'after');
            end
        
        end
        
        %% plotBlock
        function plotBlock(this)
        % Plots the block to current axis
        
            % Add last first point to last to close contour
            plot([this.keypointsX; this.keypointsX(1)],...
                [this.keypointsY; this.keypointsY(1)],'-x');
            axis equal;
            grid on
        end
    end
    
end    