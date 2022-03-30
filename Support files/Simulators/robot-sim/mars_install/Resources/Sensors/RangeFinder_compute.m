function [EXP]=RangeFinder_compute(EXP,n_robot,n_sensor)

EXP.Agent(n_robot).Sensor(n_sensor).Presence=false;

span=EXP.Agent(n_robot).Sensor(n_sensor).Angle_span;
Number_of_measures=EXP.Agent(n_robot).Sensor(n_sensor).Number_of_measures;
range=EXP.Agent(n_robot).Sensor(n_sensor).Range;
r=EXP.Robot.Diameter/2;  % raggio del robot

EXP.Agent(n_robot).Sensor(n_sensor).Measured_angle=linspace(-span/2,span/2,Number_of_measures)';
EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance=-ones(Number_of_measures,1);
G=EXP.Agent(n_robot).Sensor(n_sensor).Measured_angle;

x0=EXP.Pose(1,n_robot);
y0=EXP.Pose(2,n_robot);
a0=EXP.Pose(3,n_robot);
%-- Compute sensor output

%-- Check robots
for i=1:EXP.Robots  % ciclo sui robots
    if (i~=n_robot)   % controllo se robot � diverso da quello che sta misurando
        if norm(EXP.Geometric_center(1:2,i)-EXP.Pose(1:2,n_robot))<=range+r  % controllo se il robot � nel range del sensore
            for j=1:Number_of_measures   % ciclo sull'angolo
                x1=EXP.Geometric_center(1,i);
                y1=EXP.Geometric_center(2,i);
                xd=(x0-x1); yd=(y0-y1);
                vx=cos(a0-G(j)); vy=sin(a0-G(j));
                delta=r^2-(xd*vy-yd*vx)^2;
                gamma=(vx*xd+vy*yd);
                if (delta>=0)&&(gamma<=0)
                    distanza=-gamma-sqrt(delta);
                    %-- controllo se aggiornare la distanza
                    if (distanza<=range)&&((EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance(j)>distanza)||(EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance(j)==-1))
                        EXP.Agent(n_robot).Sensor(n_sensor).Presence=true;  % segnalo presenza di qualcosa nel campo visivo
                        EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance(j)=distanza;
                        EXP.Agent(n_robot).Sensor(n_sensor).Detected_robots = 1;
                    end
                end
            end
        end
    end
end

%-- Check obstacles
if isfield(EXP,'Map')&& isfield(EXP.Map,'Obstacle_distance')
    warning off MATLAB:nearlySingularMatrix
    for s=1:length(EXP.Map.Obstacle)  % ciclo sugli ostacoli
        if EXP.Map.Obstacle_distance(n_robot,s).Min_dist<=range % controllo che l'ostacolo sia dentro il range
            
            V=[EXP.Map.Obstacle(s).Vertex; EXP.Map.Obstacle(s).Vertex(1,1:2)];  % aggiungo il primo vertice alla fine
            
            for k=1:length(V)-1 % ciclo sui lati dell'ostacolo
                x2=V(k,1);
                x3=V(k+1,1);
                y2=V(k,2);
                y3=V(k+1,2);
                
                for j=1:Number_of_measures % ciclo sull'angolo
                    A=[cos(a0-G(j)), -(x3-x2); sin(a0-G(j)), -(y3-y2)];
                    B=[x2-x0; y2-y0];
                    X=A\B;
                    t=X(1);
                    fi=X(2);
                    
                    if (0<=t)&&(t<=range)&&(0<=fi)&&(fi<=1)
                        distanza=t;
                        
                        %-- controllo se aggiornare la distanza
                        if (distanza<=range)&&((EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance(j)>distanza)||(EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance(j)==-1))
                            EXP.Agent(n_robot).Sensor(n_sensor).Presence=true;  % segnalo presenza di qualcosa nel campo visivo
                            EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance(j)=distanza;
                            EXP.Agent(n_robot).Sensor(n_sensor).Detected_obstacles = 1;
                        end
                    end
                end
            end
        end
    end
    warning on MATLAB:nearlySingularMatrix
end

%-- Check bounds
warning off MATLAB:nearlySingularMatrix
for s=1:1  % ciclo sugli ostacoli
    %if EXP.Map.Obstacle_distance(n_robot,s).Min_dist<=range % controllo che l'ostacolo sia dentro il range
    V= EXP.Workspace;%[EXP.Map.Obstacle(s).Vertex; EXP.Map.Obstacle(s).Vertex(1,1:2)];  % aggiungo il primo vertice alla fine
    
    for k=1:length(V)-1 % ciclo sui lati dell'ostacolo
        x2=V(k,1);
        x3=V(k+1,1);
        y2=V(k,2);
        y3=V(k+1,2);
        
        for j=1:Number_of_measures % ciclo sull'angolo
            A=[cos(a0-G(j)), -(x3-x2); sin(a0-G(j)), -(y3-y2)];
            B=[x2-x0; y2-y0];
            X=A\B;
            t=X(1);
            fi=X(2);
            
            if (0<=t)&&(t<=range)&&(0<=fi)&&(fi<=1)
                distanza=t;
                
                %-- controllo se aggiornare la distanza
                if (distanza<=range)&&((EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance(j)>distanza)||(EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance(j)==-1))
                    EXP.Agent(n_robot).Sensor(n_sensor).Presence=true;  % segnalo presenza di qualcosa nel campo visivo
                    EXP.Agent(n_robot).Sensor(n_sensor).Measured_distance(j)=distanza;
                    EXP.Agent(n_robot).Sensor(n_sensor).Detected_wall = 1;
                end
            end
        end
    end
    % end
end
warning on MATLAB:nearlySingularMatrix

%------
