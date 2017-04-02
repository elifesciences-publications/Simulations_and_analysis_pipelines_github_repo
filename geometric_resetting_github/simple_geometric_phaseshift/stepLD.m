function new_phase = stepLD( theta, R, X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
xL = 0;
yL = 0;
rL = 1;

% WLOG, night limit cycle has radius R and is displaced along the x-axis X
xD = X;
yD = 0;
rD = R;

x0 = xL + rL * cos(theta*2*pi/360);
y0 = yL + rL * sin(theta*2*pi/360);
    
    % assuming instantaneous attraction to dark limit cycle, find new angle
    xjump = x0 - xD;
    yjump = y0;
    
    newTheta = theta;
    if abs(yjump) > 0
       newTheta = atand(yjump/xjump);
       xf1 = xD + rD * cos(newTheta*2*pi/360);
       yf1 = yD + rD * sin(newTheta*2*pi/360);
       xf2 = xD - rD * cos(newTheta*2*pi/360);
       yf2 = yD - rD * sin(newTheta*2*pi/360); 
       % i'm too lazy to figure out a good way to decide which branch of
       % arctangent should be used, so just test which is closest to the
       % original point.
       dist1 = (xf1 - x0)^2 + (yf1 - y0)^2;
       dist2 = (xf2 - x0)^2 + (yf2 - y0)^2;
       
       if dist1 < dist2
          xf = xf1;
          yf = yf1;
       else
          xf = xf2;
          yf = yf2;
          if newTheta > 0
             newTheta = newTheta + 180;
          else
             newTheta = newTheta - 180;
          end
          
       end

    end
new_phase = newTheta;
end

