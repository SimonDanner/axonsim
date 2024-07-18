function [ x,y ] = curvature( c, length )
    parts=499;
    s = sign(c);
    c=s*c;
    if c > 1
        c=1;
    end
    if c ~= 0
    
    r = 1/c*length/pi;
    U = 2*r*pi;
    angle = length/U*2*pi;
    dev=sin(angle/2)*r;
    x = -dev:(2*dev/parts):dev;
    y = -sqrt(r^2 - x.^2);
    y = y+r;
    end
    if c == 0
        x = -length/2:length/parts:length/2;
        y = zeros(1,parts+1);
    end
    y = s.*y;
end

