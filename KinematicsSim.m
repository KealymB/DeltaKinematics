%% Dimensions (mm)
global baseRadius; 
global BicepLength;
global ForearmLength;
global ForearmWidth;
global wristRadius;
global z_offset;
global forearmPadding;
global planeHeight;

global M1;
global M2;
global M3;
global currX;
global currY;
global currZ;
global coords;
global prevCoords;

M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0;
currZ = 286.0;
coords = [0 0 0];
prevCoords = [0 0 0];

baseRadius = 65; 
ForearmLength = 185;
BicepLength = ForearmLength*0.5;
ForearmWidth = 35;
wristRadius = 20;
z_offset = 0;
forearmPadding = 5; % used for intersection test. and will make the forarms bigger to account for the width of the forearms
planeHeight = 100;

%% testing Accuracy (Square)
grid on
M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0;
currZ = 250.0;

height = 80+25;
BicepLength = ForearmLength*0.7;

hold on

linearMove(50, -50, height, 1); %bottom left
linearMove(50, -50, height-5, 1); % drop pen
linearMove(50, 50, height-5, 1); % drop pen
linearMove(-50, 50, height-5, 1); % drop pen
linearMove(-50, -50, height-5, 1); % drop pen
linearMove(50, -50, height-5, 1); % drop pen

hold off

%% Testing Microsteps on accuracy and height
grid on
M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0; 
currZ = 250.0;

height = 80+25;
BicepLength = ForearmLength*0.7;

hold on

linearMove(0, 0, height-5, 1); %top
linearMove(50, 50, height-5, 1); %top

%% Testing Curviture
grid on
M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0; 
currZ = 250.0;

height = 80+25;
BicepLength = ForearmLength*0.7;

hold on

linearMove(0, -0, height-5, 1); %center
linearMove(50, -0, height-5, 1); %top
linearMove(-50, -0, height-5, 1); %bottom
linearMove(0, -0, height-5, 1); %bottom
linearMove(0, -50, height-5, 1); %left
linearMove(0, 50, height-5, 1); %right

hold off

%% Testing plane
planeHeight = 120;
BicepLength = ForearmLength*0.7;
planeWidth = 140;
wristZoffset = 25; % how much heigher the wrist needs to be than the drawing surface
zHeight = planeHeight + wristZoffset;
hasIntersection = 0; % used to check if there was an intersection.

testPoints = [[-planeWidth/2 planeWidth/2]; [-planeWidth/2 -planeWidth/2]; [planeWidth/2 planeWidth/2]; [planeWidth/2 -planeWidth/2]];

hold on 

for index = 1:length(testPoints)
    theta = delta_calcInverse(testPoints(index, 1), testPoints(index, 2), -zHeight);
    hasIntersection = hasIntersection + KinematicSim(theta(1), theta(2), theta(3), true);
end

if hasIntersection > 0
    disp("collision detected at height: " + planeHeight);
end


[x y] = meshgrid(-planeWidth/2:10:planeWidth/2); % Generate x and y data
z = zeros(size(x, 1))+planeHeight; % Generate z data
surf(x, y, z) % Plot the surface
hold off

%% testing Bezier curve
grid on
M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0;
currZ = 286.0;

s = [-50 0];
c1 = [-25 -50];
c2 = [25, 50];
e = [50 0];

[X, Y] = cBezier(s, c1, c2, e, 0.025);

hold on

linearMove(s(1), s(2), 180, 1); %bottom left
linearMove(s(1), s(2), 150, 1); % drop pen

for t=1:length(X)
    linearMove(X(t), Y(t), 150, 1);
end

hold off
%% testing Accuracy (Square)
grid on
M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0;
currZ = 250.0;

height = 80+25;
BicepLength = ForearmLength*0.7;

hold on

linearMove(-50, 50, height, 1); %bottom left
linearMove(-50, 40, height, 1); % drop pen

hold off
%% testing drawing SVG
grid on
M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0;
currZ = 286.0;

hold on

for i = 1: length(commands)
    if contains(commands(i), "JZ")
        part = str2double(erase(split(commands(i), " "), "!"));
        linearMove(currX, currY, part(2), 1);
    end

    if contains(commands(i), "LM")
        part = str2double(erase(split(commands(i), " "), "!"));
        linearMove(part(2), part(3), part(4), 1);
    end

    if contains(commands(i), "CB")
        part = str2double(erase(split(commands(i), " "), "!"));
        s = [part(2), part(3)];
        c1 = [part(4), part(5)];
        c2 = [part(6), part(7)];
        e = [part(8), part(9)];
        zHeight = part(10);
        
        [X, Y] = cBezier(s, c1, c2, e, 0.05);

        for t=1:length(X)
            linearMove(X(t), Y(t), zHeight, 1);
        end
    end
end

hold off

%% testing line
grid on
M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0;
currZ = 250;

hold on

linearMove(0, 0, 180, 1);   %pen move
linearMove(0, 0, 150, 1);   %pen drop
linearMove(60, 0, 150, 1);   %pen drop

hold off

%% Show arms
hold on

KinematicSim(M1, M2, M3, true);

hold off
%% testing line drawing accuracy
grid on
M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0;
currZ = 286.0;

hold on

r = 49.5;

for theta = 0:pi/4:pi
    x = r*sin(theta);
    y = r*cos(theta);
    linearMove(x, y, 180, 1);   %pen move
    linearMove(x, y, 150, 1);   %top middle
    linearMove(-x, -y, 150, 1); %top middle
    linearMove(-x, -y, 180, 1); %lift pen
end

hold off

%% testing circle drawing
grid on
M1 = 1.91986;
M2 = 1.91986;
M3 = 1.91986;
currX = 0.0;
currY = 0.0;
currZ = 286.0;


index = 0:0.01:2*pi;
hold on

for i = 0:4
    r = 50-i*15; % radius of 70mm
    x = r*sin(index);
    y = r*cos(index);
    
    linearMove(0, 50, 180, 1); %position pen top center
    
    for pos = 1:length(index)
        linearMove(x(pos), y(pos), 150, 1); % drop pen
    end
end
hold off

%% Linear move
function lm = linearMove(xt, yt, zt, stepSize)
    global M1;
    global M2;
    global M3;

    global currX;
    global currY;
    global currZ;

    global coords;
    global prevCoords;

    rad2step = 12800/(2*pi);
    step2rad = 1/rad2step; %this is the smallest "degree" the motor can turn
    
    xDist = xt - currX;
    yDist = yt - currY;
    zDist = zt - currZ;
    
    totDist = sqrt(xDist^2+yDist^2+zDist^2);
    numberOfSteps = round(totDist / stepSize);
    
    xStep = xDist / numberOfSteps;
    yStep = yDist / numberOfSteps;
    zStep = zDist / numberOfSteps;
    
    for i = 1: numberOfSteps
        xInterp = currX + i * xStep;
        yInterp = currY + i * yStep;
        zInterp = currZ + i * zStep;

        % Z calculation

        %
        % A-----B
        % |     |
        % |     |   
        % D-----C
        %, 6.0 7.0 8.0

        Aheight = 0; % location(-xMax, yMax)
        Bheight = 0 ; % location(xMax, yMax)
        Cheight = 0; % location(xMax, -yMax)
        Dheight = 0; % location(-yMax, -yMax)

        xMax = 50;
        yMax = 50; 
        
        qa = ((-yMax -yInterp)/(2*(-yMax)))*Aheight + ((yInterp - yMax)/(2*(-yMax)))*Dheight;
        qb = ((-yMax -yInterp)/(2*(-yMax)))*Bheight + ((yInterp - yMax)/(2*(-yMax)))*Cheight;

        %zHeight = (-xMax -xInterp)/(2*(xMax))*qa + (xInterp - xMax)/(2*(xMax))*qb;

        radius = sqrt(xInterp^2 + yInterp^2);

        addHeight = -1*abs((70.71-radius)/70.71);

        zHeight = 0;
    
        theta = delta_calcInverse(xInterp, yInterp, -(zHeight+zt));
    
        error = mod(theta, step2rad);
        theta = theta - error; %this removes the left over angle? 
    
        coords = KinematicSim(theta(1), theta(2), theta(3), false);

        if prevCoords(1) ~= 0 && prevCoords(2) ~= 0 && prevCoords(3) ~= 0
            if zInterp <= 150
                plot3([prevCoords(1) coords(1)], [prevCoords(2) coords(2)], [prevCoords(3) coords(3)], "black");
            end
        end
%         plot3(coords(1), coords(2), coords(3), "o");
    
        M1 = theta(1);
        M2 = theta(2);
        M3 = theta(3);
        prevCoords = coords;
    end

    currX = coords(1);
    currY = coords(2);
    currZ = coords(3);
end
%% Forward Kinematics Simulation
function coords = KinematicSim(th1, th2, th3, plotArms)
    global baseRadius; 
    global BicepLength;
    global ForearmLength;
    global ForearmWidth;
    global wristRadius;
    global z_offset;
    global planeHeight;
    
    % utils for base
    angles = [0 120 240]*(pi/180);
    x_actuator = baseRadius*sin(angles);
    y_actuator = baseRadius*cos(angles);
    z_actuator = [0 0 0];
    
    %  B1
    y1 = y_actuator(1)-BicepLength*cos(-th1);
    z1 = z_actuator(1)-BicepLength*sin(-th1);
    
    x_b1 = [x_actuator(1) 0];
    y_b1 = [y_actuator(1) y1];
    z_b1 = [z_actuator(1) z1];
    
    % B2
    x2 = x_actuator(2)-(sqrt(3)/2)*BicepLength*cos(-th2);
    y2 = y_actuator(2)+0.5*BicepLength*cos(-th2);
    z2 = z_actuator(2)-BicepLength*sin(-th2);
    
    x_b2 = [x_actuator(2) x2];
    y_b2 = [y_actuator(2) y2];
    z_b2 = [z_actuator(2) z2];
    
    % B3
    x3 = x_actuator(3)+(sqrt(3)/2)*BicepLength*cos(-th3);
    y3 = y_actuator(3)+0.5*BicepLength*cos(-th3);
    z3 = z_actuator(3)-BicepLength*sin(-th3);
    
    x_b3 = [x_actuator(3) x3];
    y_b3 = [y_actuator(3) y3];
    z_b3 = [z_actuator(3) z3];
    

    if plotArms
        plot3(x_b1,y_b1,z_b1, "Color","r");
        plot3(x_b2,y_b2,z_b2, "Color","r");
        plot3(x_b3,y_b3,z_b3, "Color","r");
    end
    
    % forearms
    y1 = y1-wristRadius;
    z1 = BicepLength*sin(th1);
    % plot3([0 0],[y_actuator(1) y1],[z_actuator(1) z1]);
    % Sphere(0, y1, z1, ForearmLength);
    
    y2 = y2+wristRadius*sin(pi/6);
    x2 = -y2*tan(pi/3);
    z2 = BicepLength*sin(th2);
    % plot3([x_actuator(2) x2],[y_actuator(2) y2],[z_actuator(2) z2]);
    % Sphere(x2, y2, z2, ForearmLength);
    
    y3 = y3+wristRadius*sin(pi/6);
    x3 = y3*tan(pi/3);
    z3 = BicepLength*sin(th3);
    % Sphere(x3, x3, z3, ForearmLength);
    
    % intersection testing
    points = intersectionPoints([0 y1 z1], [x2 y2 z2], [x3 y3 z3], ForearmLength);
    x0 = points(4);
    y0 = points(5);
    z0 = points(6);
    
    % F1
    x_f1 = [0 x0];
    y_f1 = [y_actuator(1)-BicepLength*cos(-th1) y0+wristRadius];
    z_f1 = [z_actuator(1)-BicepLength*sin(-th1) z0];
%     intersecting = LinePlaneIntersection([0 y_actuator(1)-BicepLength*cos(-th1) z_actuator(1)-BicepLength*sin(-th1)], [x0 y0+wristRadius z0], [0 0 planeHeight+10], [0 0 1]);
    
    % F2
    x_f2 = [x_actuator(2)-(sqrt(3)/2)*BicepLength*cos(-th2) x0+wristRadius*cos(pi/6)];
    y_f2 = [y_actuator(2)+0.5*BicepLength*cos(-th2) y0-wristRadius*sin(pi/6)];
    z_f2 = [z_actuator(2)-BicepLength*sin(-th2) z0];
%     intersecting = intersecting + LinePlaneIntersection([x_actuator(2)-(sqrt(3)/2)*BicepLength*cos(-th2) y_actuator(2)+0.5*BicepLength*cos(-th2) z_actuator(2)-BicepLength*sin(-th2)], [x0+wristRadius*cos(pi/6) y0-wristRadius*sin(pi/6) z0], [0 0 planeHeight+10], [0 0 1]);
    
    % F3
    x_f3 = [x_actuator(3)+(sqrt(3)/2)*BicepLength*cos(-th3) x0-wristRadius*cos(pi/6)];
    y_f3 = [y_actuator(3)+0.5*BicepLength*cos(-th3) y0-wristRadius*sin(pi/6)];
    z_f3 = [z_actuator(3)-BicepLength*sin(-th3) z0];
%     intersecting = intersecting + LinePlaneIntersection([x_actuator(3)+(sqrt(3)/2)*BicepLength*cos(-th2) y_actuator(3)+0.5*BicepLength*cos(-th2) z_actuator(3)-BicepLength*sin(-th2)], [x0-wristRadius*cos(pi/6) y0-wristRadius*sin(pi/6) z0], [0 0 planeHeight+10], [0 0 1]);
    
    if plotArms
        plot3(x_f1,y_f1,z_f1, "Color","b");
        plot3(x_f2,y_f2,z_f2, "Color","b");
        plot3(x_f3,y_f3,z_f3, "Color","b");
    end

    if plotArms
        plot3(x0,y0,z0, "o"); % point from reverse kinematics
    end
%     hasIntersection = intersecting;
    coords = [x0, y0, z0];
end
%% Inverse Kinematic Functions
function AngleYZ = delta_calcAngleYZ(x_, y_, z_)
    global baseRadius; 
    global wristRadius;
    global BicepLength;
    global ForearmLength;

    baseR = baseRadius + 0;
    wristR = wristRadius + 0;
    bicepL = BicepLength + 0;
    ForearmL = ForearmLength +0;
    z_ = z_ + 0;

    arm_end_y = y_+wristRadius; % add wrist radius when designing more of the arm aka: arm_end_x = x_ +wristRadius;
    l2_XZ = sqrt(ForearmLength^2 - x_^2); %The length of link 2 when projected onto the YZ plane

    l2_angle = asin(y_/ForearmLength); % can use this to avoid too angled position due to ball joints
    
    ext = sqrt(z_^2 + (baseR - arm_end_y)^2); %Extension of the arm from the centre of the servo rotation to the end ball joint of link2
    
    phi = acos((BicepLength^2 + ext^2 - l2_XZ^2) / (2 * BicepLength * ext)); %Cosine rule that calculates the angle between the ext line and L1
    omega = -atan2(z_, baseR - arm_end_y); %Calculates the angle between horizontal (X) the ext line with respect to its quadrant
    theta = phi + omega; %Theta is the angle between horizontal (X) and L1
    AngleYZ = theta;
end

function angles = delta_calcInverse(x_, y_, z_)
    theta1 = 0;
    theta2 = 0;
    theta3 = 0;

    status = delta_calcAngleYZ(x_, y_, z_);
    if (status ~= -1)
        theta1 = status;
    end
    status = delta_calcAngleYZ(x_*cos(2*pi/3) - y_*sin(2*pi/3), y_*cos(2*pi/3)+x_*sin(2*pi/3), z_);
    if (status ~= -1)
        theta2 = status;
    end
        status = delta_calcAngleYZ(x_*cos(2*pi/3) + y_*sin(2*pi/3), y_*cos(2*pi/3)-x_*sin(2*pi/3), z_);
    if (status ~= -1)
        theta3 = status;
    end
    angles = [theta1 theta2 theta3];
end

function distance = distanceCalc(x,y,z)
    distance = sqrt((x(1) - x(2))^2+(y(1)- y(2))^2+(z(1)- z(2))^2);
end

function intersection = intersectionPoints(p1, p2, p3, r)
    p21 = p2-p1;
    p31 = p3-p1;
    c = cross(p21,p31);
    c2 = sum(c.^2);
    u1 = cross(((sum(p21.^2)+r^2-r^2)*p31 - ...
             (sum(p31.^2)+r^2-r^2)*p21)/2,c)/c2;
    v = sqrt(r^2-sum(u1.^2))*c/sqrt(c2);
    i1 = p1+u1+v;
    i2 = p1+u1-v;
    intersection = [i1, i2];
end

%% Sim functions

function Sphere(x,y,z, r)
    [X,Y,Z] = sphere;
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    
    s = surf(X2+x,Y2+y,Z2+z,"FaceAlpha",0.2);
    s.EdgeColor = "none";

    % plot center point
    plot3(x,y,z, "o");
end

function isIntersecting = LinePlaneIntersection(L_start, L_end, P_point, P_normal, forearmPadding)
    global forearmPadding;
    L_normal = L_end-L_start;
    P_point(3) = P_point(3)-10;
    [I,rc] = line_plane_intersection(L_normal, L_start, P_normal, P_point);

    if abs(I(2)) >= 70+forearmPadding
        isIntersecting = false;
    elseif abs(I(1)) >= 70+forearmPadding
        isIntersecting = false;
    else
        isIntersecting = (rc > 0);
    end
end


%% Bezier functions
function [X, Y] = cBezier(p0, p1, p2, p3, interpDist)
    START = p0;
    C1 = p1;
    C2 = p2;
    END = p3;
        
    x = [];
    y = [];

    for T=0:interpDist:1
        %quad bezier
        xa = lerp( START(1) , C1(1) , T );
        ya = lerp( START(2) , C1(2) , T );
        xb = lerp( C1(1) , C2(1) , T );
        yb = lerp( C1(2) , C2(2) , T );
        xc = lerp( C2(1) , END(1) , T );
        yc = lerp( C2(2) , END(2) , T );
        
        %cubic bezier
        xm = lerp( xa , xb , T );
        ym = lerp( ya , yb , T );
        xn = lerp( xb , xc , T );
        yn = lerp( yb , yc , T );

        %path coords
        x(end+1) = lerp(xm, xn, T);
        y(end+1) = lerp(ym, yn, T);
    end

    X = x;
    Y = y;
end

function linearInterp = lerp(v1, v2, t)
    linearInterp = v1 + (v2-v1)*t;
end