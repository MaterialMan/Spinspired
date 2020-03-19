%% Two-dimensional Izhikevich model 
clear

% 1 ) initialise parameters  
dt = 0.5;

a = 0.02;
b = 0.2;
c = -65;
d = 2;%8;

v_t = 30; %threshold mV
T_span = 100; %t/ms
T = ceil(T_span/dt) ;
I_input = 20; %inject current nA
I_t = round(rand(T,1)); % input timing

% 2 ) Reserve memory
v = zeros(T,1);
u = zeros(T,1);
Iapp = zeros(T,1);
v(1) = -70; % resting potential
u(1) = -14; % steady state

%3 ) for loop
for t = 1:T-1
    
    %3.1) get input
    if I_t(t)%t*dt >10 && t*dt <90
        Iapp(t) = I_input; % synaptic current injected into neuron
    else
        Iapp(t) = 0;
    end
    
    if v(t)<v_t
        %3.2) update ODE
        dv = (0.04*v(t)+5)*v(t)+140-u(t);
        v(t+1) = v(t) + (dv+Iapp(t))*dt; % membrane potential of neuron
        
        du = a*(b*v(t)-u(t));
        u(t+1) = u(t) + dt*du ; % recovery variable
    else
        %3 . 3 ) spike!
        v(t) = v_t;
        v(t+1) = c;
        u(t+1) = u(t)+d ;
    end
end

% 4 ) plot voltage trace
figure
subplot(1,2,1)
plot((0:T-1)*dt,v,'b') ;
xlabel('Time [ ms ]') ;
ylabel('Membrane voltage [mV]') ;

subplot(1,2,2)
plot((0:T-1)*dt,Iapp,'b') ;
xlabel('Time [ ms ]') ;
ylabel('Membrane voltage [mV]') ;
