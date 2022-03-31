
function [x,y, num_points, area_ratio] = getPolyShape(shape, num_sides, rotate_angle, ref_point)

start_pos = 0;

switch(shape)
    case 'custom'
        % create array containing 361 equally spaced points btw 0 and 2*pi
        degrees = linspace(0, 2*pi, 361);
        % store in the array 'c' cosine of all the values in the array 'degrees'
        c = cos(degrees);
        % store in the array 's' sine of all the values in the array 'degrees'
        s = sin(degrees);
        % calculate appropriate step size for plotting a hexagon
        step = 360/num_sides;
        %plot the polygon
        x = (c(1:step:361)+1)/2;
        y = (s(1:step:361)+1)/2;
               
%     case 'donut'
%         t = 0.05:0.25:2*pi;
%         x1 = cos(t);
%         y1 = sin(t);
%         x2 = 0.5*cos(t);
%         y2 = 0.5*sin(t);
%         pgon = polyshape({(x1+1)/2,(x2+1)/2},{(y1+1)/2,(y2+1)/2});
%         %pgon = polyshape({x1,x2},{y1,y2});
%         x = pgon.Vertices(:,1)';
%         y = pgon.Vertices(:,2)';
%     
%         area_ratio = polyarea((x1+1)/2,(y1+1)/2);
%         
%     case '4-square'
%         % create 4 separate squares 
%         pgon(1) = nsidedpoly(4,'Center',[0.25 0.75],'SideLength',0.25);
%         pgon(2) = nsidedpoly(4,'Center',[0.25 0.25],'SideLength',0.25);
%         pgon(3) = nsidedpoly(4,'Center',[0.75 0.25],'SideLength',0.25);
%         pgon(4) = nsidedpoly(4,'Center',[0.75 0.75],'SideLength',0.25);
%         
%         x = []; y =[];
%         for i = 1:4
%             x = [x pgon(i).Vertices(:,1)'];
%             y = [y pgon(i).Vertices(:,2)'];
%         end
%         
    case 'arrowhead'
        x = [0.1 0.4 0.4 0.9 0.9 0.4 0.4];
        y = [0.5 0.85 0.7 0.7 0.3 0.3 0.15];
    case 'square'
        x = [start_pos 0 1 1];
        y = [start_pos 1 1 0];
    case 'rectangle'
        x = [start_pos 0    1    1];
        y = [start_pos 0.25  0.25  0];
    case 'triangle'
        x = [start_pos 0 1 0.5];
        y = [start_pos 0 0 1];
    case 'trapezoid'
        x = [start_pos 0.2 0.8 1];
        y = [start_pos 0.5 0.5 0];
    case 'parallelogram'
        x = [start_pos 0.2 1 0.8];
        y = [start_pos 0.5 0.5 0];
%     case 'pentogram'
%         x = [start_pos 0.2 1 0.8];
%         y = [start_pos 0.5 0.5 0];
    case 'circle'
        num_points = 37;
        m = [1                  0.500000000000000
            0.992403876506104	0.586824088833465
            0.969846310392954	0.671010071662834
            0.933012701892219	0.750000000000000
            0.883022221559489	0.821393804843270
            0.821393804843270	0.883022221559489
            0.750000000000000	0.933012701892219
            0.671010071662834	0.969846310392954
            0.586824088833465	0.992403876506104
            0.500000000000000	1
            0.413175911166535	0.992403876506104
            0.328989928337166	0.969846310392954
            0.250000000000000	0.933012701892219
            0.178606195156730	0.883022221559489
            0.116977778440511	0.821393804843270
            0.0669872981077807	0.750000000000000
            0.0301536896070458	0.671010071662834
            0.00759612349389599	0.586824088833465
            0                   0.500000000000000
            0.00759612349389593	0.413175911166535
            0.0301536896070458	0.328989928337166
            0.0669872981077806	0.250000000000000
            0.116977778440511	0.178606195156730
            0.178606195156730	0.116977778440511
            0.250000000000000	0.0669872981077806
            0.328989928337166	0.0301536896070458
            0.413175911166535	0.00759612349389599
            0.500000000000000	0
            0.586824088833465	0.00759612349389593
            0.671010071662834	0.0301536896070457
            0.750000000000000	0.0669872981077807
            0.821393804843270	0.116977778440511
            0.883022221559489	0.178606195156730
            0.933012701892219	0.250000000000000
            0.969846310392954	0.328989928337166
            0.992403876506104	0.413175911166535
            1	0.500000000000000];
        x = m(:,1)';
        y = m(:,2)';
end


% rotate shape
if nargin > 2
    % create poly struct
    pgon = polyshape(x,y);
    rot_pgon = rotate(pgon,rotate_angle,ref_point);
    
    x = rot_pgon.Vertices(:,1)';
    y = rot_pgon.Vertices(:,2)';
    
    % plot shape rotation
    %subplot(1,2,1)
    %plot(pgon)
    %subplot(1,2,2)
    %plot(rot_pgon)
    %drawnow
end
% 
num_points = length(x); 

if max(x) > 1 || min(x) < 0 || max(y) > 1 || min(y) < 0
    error('Vertices out of bounds')
end

if ~exist('area_ratio')
area_ratio = polyarea(x,y);
end