function draw_robot(x,y,heading_angle)
persistent nverts npts verts linewidth color hrobot sensorpt
global new_verts fverts xR yR xL yL

%persistent variables are local to the function but kept in memory for the
%next function call

if(nargin < 3)
    % initialize
    [verts, linewidth, color] = robot_coords;
    fverts=verts;
    nverts = size(verts,1);
    npts = size(verts, 2);
    figure(1); hold on;
    % Show the robot's current location
    for i=1:nverts
        hrobot(i) = plot(verts(i,1:(npts/2)), verts(i,(npts/2+1):npts),'color',color(i,:),'linewidth',linewidth(i));
    end
    
else
    % redraw the robot in the new configuration
    tmp_verts = verts;
    
    % adjust heading - rotate the entire robot
    cosa = cos(heading_angle);
    sina = sin(heading_angle);
    new_verts(:,1:(npts/2)) =  x + cosa * tmp_verts(:,1:(npts/2)) - sina * tmp_verts(:,(npts/2+1):npts);
    new_verts(:,(npts/2+1):npts) =  y + sina * tmp_verts(:,1:(npts/2)) + cosa * tmp_verts(:,(npts/2+1):npts);
    
    figure(1); hold on
    for i=1:nverts
        set(hrobot(i),'xdata',new_verts(i,1:(npts/2)),'ydata',new_verts(i,(npts/2+1):npts));
    end
    
end
%--------------------------------------------------------------------------
% end of draw_robot
%--------------------------------------------------------------------------