function [wheelimg rangle gangle bangle] = generate_colorwheel(rad)


wheelimg = zeros(2*rad+1,2*rad+1,3);
newwheel = zeros(size(wheelimg,1)+2*round(rad/5),size(wheelimg,2)+2*round(rad/5),3);
rise = linspace(0,1,61);
fall = linspace(1,0,61);

rangle = 0.85*[ones(1,59),fall,zeros(1,119),rise,ones(1,60)]+0.1;
gangle = 0.85*[rise(2:end),ones(1,59),ones(1,60),fall,zeros(1,120)]+0.1;
bangle = 0.85*[zeros(1,119),rise,ones(1,59),ones(1,60),fall]+0.1;

rangle = [rangle, rangle];
gangle = [gangle, gangle];
bangle = [bangle, bangle];


for m = 1:2*rad+1
    for n = 1:2*rad+1
        x = (n-(rad+1));
        y = -(m-(rad+1));
        unitmag = (1/rad)*(x^2 + y^2)^0.5;
        angle = round(atand(y/x));
        if x < 0
            angle = 180 + angle; 
        elseif y < 0 && x >= 0
            angle = 360 + angle;
        end
        if angle == 0
            angle = 360;
        end
        if unitmag <= 1 && unitmag >=0.5 && unitmag ~= 0
            wheelimg(m,n,1) = rangle(1,angle);%unitmag*rangle(1,angle)+(1-unitmag)*0.25;
            wheelimg(m,n,2) = gangle(1,angle);%unitmag*gangle(1,angle)+(1-unitmag)*0.25;
            wheelimg(m,n,3) = bangle(1,angle);%unitmag*bangle(1,angle)+(1-unitmag)*0.25;
%         elseif unitmag == 0
%             wheelimg(m,n,:) = 0.25*ones(1,1,3);
        end
    end
end

newwheel(round(rad/5)+1:round(rad/5)+size(wheelimg,1),round(rad/5)+1:round(rad/5)+size(wheelimg,2),:) = wheelimg;

newwheel(rad+round(rad/5)+1,:,:) = 1;
newwheel(:,rad+round(rad/5)+1,:) = 1;

wheelimg = newwheel;