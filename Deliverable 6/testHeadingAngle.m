clear; clc

heading_angle = 3*pi/2;
x = 0; y = 0;

vL = [0 50 10 10 50];
vR = [0 10 50 50 10];

for k = 1:length(vL)    
    vel_left = vL(k);
    vel_right = vR(k);
    
    if vel_left > vel_right
        heading_angle = heading_angle - pi/2 + atan(10/(vel_left-vel_right));
    elseif vel_right > vel_left
        heading_angle = heading_angle + atan( (vel_right-vel_left)/10 );
    end
    
    figure(1)
    draw_robot(x,y, heading_angle);
    keyboard
end