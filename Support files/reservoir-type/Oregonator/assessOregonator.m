%% assessOregonator.m
% Function to collect reservoir states for oregonator reservoir.
%
% notes: 
% - reactor(s) and inputs are built and assigned before running through samples
% - inputs have a time-period. After n time-steps the next task input
% sample is given
% - States are currently read as convolved summation of pixels in vesicle

function[final_states,individual] = assessOregonator(individual,input_sequence,config,target_output)

%if single input entry, add previous state
if size(input_sequence,1) == 1
    input_sequence = [zeros(size(input_sequence)); input_sequence];
end

max_input_length = 0;

% pre-allocate state matrices
for i= 1:config.num_reservoirs
    if size(input_sequence,1) == 2
        states{i} = individual.last_state{i};
    else
        states{i} = zeros(size(input_sequence,1),individual.nodes(i));
    end
    x{i} = zeros(size(input_sequence,1),individual.nodes(i));
end

% pre-assign anything that can be calculated before running the reservoir
for i= 1:config.num_reservoirs
    [reactor{i},pos_x,pos_y]= createReactor(individual.bitmatrix{i}, individual.phimatrix{i}, individual, config);
    reactor_state{i} = zeros(size(reactor{i},1),size(reactor{i},2),4);%zeros(individual.height, individual.width, 4); %initReactor(reactor{i}, config);
    reactor_state{i}(:,:,config.p_idx) = reactor{i}(:,:,2);
    reactor_state{i}(:,:,config.b_idx)= reactor{i}(:,:,1);

    id(i,:) = pos_x;%config.vesicle_radius:config.vesicle_radius*2:2*config.vesicle_radius*individual.vesicle_grid_size(i);
end


% time multiplex -
% extra_time = 20;
% input_sequence = [zeros(extra_time,config.task_num_inputs); input_sequence];
input{i} = [[zeros(1,config.task_num_inputs); input_sequence] repmat(individual.bias_node,size(input_sequence,1)+1,1)]*(individual.input_weights{i}*individual.input_scaling(i))';
input_mul{i} = zeros((size(input_sequence,1)+1)*individual.time_period(i),size(input{i},2));
if individual.time_period(i) > 1
    for p = 0:individual.input_length(i)-1
        input_mul{i}(mod(1:size(input_mul{i},1),individual.time_period(i)) == p,:) = input{i};
    end
else
    input_mul{i} = input{i};
end

if size(input_mul{i},1) > max_input_length
    max_input_length = size(input_mul{i},1)-individual.time_period(i);
end

for n = 2:max_input_length
    in_grid = reshape(input_mul{i}(n,:),individual.vesicle_grid_size(i),individual.vesicle_grid_size(i));
    
    in_grid = in_grid.*(individual.bitmatrix{i} ==1);
    
    empty_grid = zeros(size(reactor_state{i},1));
    empty_grid(id(i,:),id(i,:)) = in_grid;
    input_grid(:,:,n) = adjustInputShape(empty_grid,individual.input_widths(i));
end
% pre-assign
%maxPhiValue = max(individual.phimatrix, [], 'all');
colormap(config.figure_array(3),'hot')

% Calculate reservoir states - general state equation for multi-reservoir system: x(n) = f(Win*u(n) + S)
for n = 2:max_input_length
    
    for i= 1:config.num_reservoirs % cycle through sub-reservoirs
                   
        % supply source inputs
        %if mod(n,individual.time_period(i)) == 0
        %sourcematrix = reactor_state{i}(id(i,:),id(i,:),config.u_idx);
        %sourcematrix(individual.sourcematrix{i} == 1) = 1;
        %reactor_state{i}(id(i,:),id(i,:),config.u_idx) = sourcematrix;%reactor_state{i}(id(i,:),id(i,:),config.u_idx)
        %reactor_state{i}(:,:,config.u_idx) = adjustInputShape(reactor_state{i}(:,:,config.u_idx),individual.input_widths(i));
        % end
        
        % change light intensity
        reactor_state{i}(:,:,config.u_idx) = reactor_state{i}(:,:,config.u_idx) + input_grid(:,:,n);

        % run sim step
        reactor_state{i} = oreg_step(reactor_state{i},size(reactor_state{i},1),size(reactor_state{i},2), config);
               
        % convolve states
        conv_state = reactor_state{i}(:,:,1);%sum(reactor_state{i},3);
        if individual.kernel_stride(i) ~= 1
            conv_state = conv2(padarray(conv_state, [individual.pad_size(i) individual.pad_size]), individual.kernel{i}, 'valid');
        end
        conv_state = conv_state(1:individual.kernel_stride(i):end, 1:individual.kernel_stride(i):end);
        
        %conv_state(isnan(conv_state)) = 0; 
        has_nan_states = nnz(isnan(conv_state)) > 2;
        conv_state(conv_state > 1) = 1;
        conv_state(conv_state < 0) = 0;
        
        if mod(n, config.stride) == 0 && config.plot_states       
            updateDisplay(reactor_state{i},size(reactor_state{i},1),size(reactor_state{i},2), config);
            title(strcat('Time:',' ',num2str(n/config.stride)));%individual.time_period(i))))
            set(0,'currentFigure',config.figure_array(3))
            %subplot(1,2,1)
            imagesc(conv_state)
            colorbar
            caxis([0 1])
            
%             subplot(1,2,2)
%             imagesc(individual.sourcematrix{i})
%             colorbar
%             caxis([0 1])
%             
%             subplot(2,2,3)
%             imagesc(reactor_state{i}(:,:,1))
%             colorbar
%             caxis([0 1])
            
            drawnow
        end
        
        % check to if simulation fails
        if has_nan_states% if failed
            states{i} = zeros(max_input_length,individual.nodes(i));
            flag =1;
            break;
        else
            flag =0;
            states{i}(n,:) = conv_state(:);
        end
        
    end
    
    if(flag==1)
        break
    end
end

%need to check! deplex to get states
for i= 1:config.num_reservoirs
    if individual.time_period(i) > 1
        states{i} = states{i}(mod(1:size(states{i},1),individual.time_period(i)) == 0,:);
    end
end

% Add leak states, if used
if config.leak_on
    states = getLeakStates(states,individual,input_sequence,config);
end

% Concat all states for output weights
final_states = [];
for i= 1:config.num_reservoirs
    final_states = [final_states states{i}];
    
    %assign last state variable
    individual.last_state{i} = states{i}(end,:);
end



% Concat input states
if config.add_input_states == 1
    final_states = [final_states input_sequence];
end

% Remove washout and output final states
if size(input_sequence,1) == 2
    final_states = final_states(end,:); % remove washout
else
    final_states = final_states(config.wash_out+1:end,:); % remove washout
end

% can remove: used as quick method to view internal reservoir activity
%if config.plot_states 
set(0,'currentFigure',config.figure_array(2))
plot(final_states)
drawnow
%end

end

%% create reactor from matrix
function [reactor,pos_x,pos_y] = createReactor(bitMatrix, phiMatrix, individual,config)

% takes the bit matrix and uses it to define the corresponding connected vesicle reactor
% config needs size of vesicle
margin_w = 3;
margin_h = 3;

reactor = zeros(individual.height+margin_h, individual.width+margin_w, 2);

vesicle_w = config.vesicle_radius * 2;
vesicle_h = vesicle_w;

% calculate inter-vesicle distance (IVD)

crossover_dist = config.crossover_dist;
IVD = vesicle_w - crossover_dist;

% calculate size of vesicle `domain' (i.e the area in the reactor within which we will place the vesicles)
m_height = size(bitMatrix,1);
m_width  = size(bitMatrix,2);

d_width  = m_width * IVD + (margin_w * 2);
d_height = m_height * IVD + (margin_h * 2);

if d_width > size(reactor,1)
    error(strcat('Vesicle domain: ', string(d_width), ' is too wide for reactor: ', string(size(reactor,1))));
end

if d_height > size(reactor,2)
    error('Vesicle domain is too high for reactor');
end

% distribute vesicle centres across vesicle domain, at IVD intervals

% first vesicle centre: (x, y)
v_centre_x = margin_w + config.vesicle_radius;
v_centre_y = margin_h + config.vesicle_radius;

% use bitmatrix to draw vesicles around desired vesicle centres
for y = 1:m_height
    v_centre_x = margin_w + config.vesicle_radius;
    pos_y(y) = v_centre_y;
    for x = 1:m_width
        pos_x(x) = v_centre_x;        
        if bitMatrix(y, x) == 1
            % draw vesicle
            reactor = drawVesicle(reactor, v_centre_x, v_centre_y, config.vesicle_radius);
        end
        v_centre_x = v_centre_x + IVD;
        
    end
    v_centre_y = v_centre_y + IVD;
    
end

%imagesc(reactor(:,:,1))

% remove vesicle boundary from adjacent vesicles
v_centre_y = margin_h + config.vesicle_radius;
for y = 1:m_height
    v_centre_x = margin_w + config.vesicle_radius;
    for x = 1:m_width
        left = false;
        top = false;
        
        % check if vesicle to the left or top
        if x > 1
            % check left
            if bitMatrix(y, x - 1) == 1 && bitMatrix(y, x) == 1
                left = true;
            end
        end
        if y > 1
            if bitMatrix(y - 1, x) == 1 && bitMatrix(y, x) == 1
                top = true;
            end
        end
        % remove overlapping boundary if present
        phi = 0;
        if bitMatrix(y, x) == 1
            phi = phiMatrix(y, x);
        end
        
        if left
            reactor = makeConnection(reactor, v_centre_x, v_centre_y, v_centre_x - IVD, v_centre_y, config.vesicle_radius);
            reactor = makeConnection(reactor, v_centre_x - IVD, v_centre_y, v_centre_x, v_centre_y, config.vesicle_radius);
            
            reactor = setExcitability(reactor, phi, v_centre_x, v_centre_y, v_centre_x - IVD, v_centre_y, config.vesicle_radius);
        end
        
        if top
            reactor = makeConnection(reactor, v_centre_x, v_centre_y, v_centre_x, v_centre_y - IVD, config.vesicle_radius);
            reactor = makeConnection(reactor, v_centre_x, v_centre_y - IVD, v_centre_x, v_centre_y, config.vesicle_radius);
            
            reactor = setExcitability(reactor, phi, v_centre_x, v_centre_y, v_centre_x, v_centre_y - IVD, config.vesicle_radius);
        end
        
        if ~left && ~top
            reactor = setExcitability(reactor, phi, v_centre_x, v_centre_y, 0, 0, config.vesicle_radius);
        end
        
        
        v_centre_x = v_centre_x + IVD;
    end
    v_centre_y = v_centre_y + IVD;
end

%imagesc(reactor(:,:,1))
% done.

% later: set phi for each vesicle

end

%% set the phi for the given vesicle
function reactor = setExcitability(reactor, phi, x0, y0, x1, y1, radius)

% find all pixels within the vesicle
% if there are connected vesicles,
%   find the locus between (x0, y0) and (x1, y1)
%   this locus forms the boundary between the two vesicles
%   remove pixels from array that are closer to (x1, y1) than (x0, y0)
% set phi for all pixels in array

hasConnection = false;
if x1 > 0 && y1 > 0
    hasConnection = true;
end

% draw box around the vesicle:
for y = (y0 - radius):(y0 + radius)
    for x = (x0 - radius):(x0 + radius)
        if euclDist(x, y, x0, y0) <= radius
            
            if hasConnection
                if euclDist(x, y, x0, y0) < euclDist(x, y, x1, y1)
                    reactor(y, x, 2) = phi;
                end
            else
                reactor(y, x, 2) = phi;
            end
            
        end
    end
end

end

function dist = euclDist(x_0,y_0,x_1,y_1)

dx = x_0 - x_1;
dy = y_0 - y_1;

dist = sqrt(dx * dx + dy * dy);

end

%% form connection between vesicles
function reactor = makeConnection(reactor, x0, y0, x1, y1, radius)

% remove all points on the circle which are < radius from the other centre-point

f = 1 - radius;
ddF_x = 0;
ddF_y = -2 * radius;
x = 0;
y = radius;

if euclDist(x0, y0 + radius, x1, y1) < radius
    reactor(y0 + radius, x0, 1) = 0.0;
end

if euclDist(x0, y0 - radius, x1, y1) < radius
    reactor(y0 - radius, x0, 1) = 0.0;
end

if euclDist(x0 + radius, y0, x1, y1) < radius
    reactor(y0, x0 + radius, 1) = 0.0;
end

if euclDist(x0 - radius, y0, x1, y1) < radius
    reactor(y0, x0 - radius, 1) = 0.0;
end

while (x < y)
    if f >= 0
        y     = y - 1;
        ddF_y = ddF_y + 2;
        f     = f + ddF_y;
    end
    
    x     = x + 1;
    ddF_x = ddF_x + 2;
    f     = f + ddF_x + 1;
    
    if euclDist(x0 + x, y0 + y, x1, y1) < radius
        reactor(y0 + y, x0 + x, 1) = 0.0;
    end
    
    if euclDist(x0 - x, y0 + y, x1, y1) < radius
        reactor(y0 + y, x0 - x, 1) = 0.0;
    end
    
    if euclDist(x0 + x, y0 - y, x1, y1) < radius
        reactor(y0 - y, x0 + x, 1) = 0.0;
    end
    
    if euclDist(x0 - x, y0 - y, x1, y1) < radius
        reactor(y0 - y, x0 - x, 1) = 0.0;
    end
    
    if euclDist(x0 + y, y0 + x, x1, y1) < radius
        reactor(y0 + x, x0 + y, 1) = 0.0;
    end
    
    if euclDist(x0 - y, y0 + x, x1, y1) < radius
        reactor(y0 + x, x0 - y, 1) = 0.0;
    end
    
    if euclDist(x0 + y, y0 - x, x1, y1) < radius
        reactor(y0 - x, x0 + y, 1) = 0.0;
    end
    
    if euclDist(x0 - y, y0 - x, x1, y1) < radius
        reactor(y0 - x, x0 - y, 1) = 0.0;
    end
    
end

end

%% draw a vesicle on the reactor
function reactor = drawVesicle(reactor, x0, y0, radius)

f = 1 - radius;
ddF_x = 0;
ddF_y = -2 * radius;
x = 0;
y = radius;

reactor(y0 + radius, x0, 1) = 1.0;
reactor(y0 - radius, x0, 1) = 1.0;
reactor(y0, x0 + radius, 1) = 1.0;
reactor(y0, x0 - radius, 1) = 1.0;


while (x < y)
    if f >= 0
        y     = y - 1;
        ddF_y = ddF_y + 2;
        f     = f + ddF_y;
    end
    
    x     = x + 1;
    ddF_x = ddF_x + 2;
    f     = f + ddF_x + 1;
    
    reactor(y0 + y, x0 + x, 1) = 1.0;
    reactor(y0 + y, x0 - x, 1) = 1.0;
    reactor(y0 - y, x0 + x, 1) = 1.0;
    reactor(y0 - y, x0 - x, 1) = 1.0;
    
    reactor(y0 + x, x0 + y, 1) = 1.0;
    reactor(y0 + x, x0 - y, 1) = 1.0;
    reactor(y0 - x, x0 + y, 1) = 1.0;
    reactor(y0 - x, x0 - y, 1) = 1.0;  
end

end

%% apply input?
function state = stimulateInput(state, amount, config)
    state(config.input_y, config.input_x, config.u_idx) = amount;
end

%% initialise reactor variables
function state = initReactor(reactor, config)

state = zeros(config.height, config.width, 4);

%     state(round(config.height / 2), round(config.width / 2), config.u_idx) = 1.0;
state(config.input_y, config.input_x, config.u_idx) = 1.0;
state(:,:,config.p_idx) = reactor(:,:,2); %config.phi_active;
state(:,:,config.b_idx) = reactor(:,:,1);

end

% udate tank view
function img = updateDisplay(state,height, width, config)

    % take state and produce image for display
    % display image in window?
  
    set(0,'currentFigure',config.figure_array(1))
   % subplot(1,2,1)
    R = state(:,:,config.u_idx);
    G = state(:,:,config.p_idx);
    B = state(:,:,config.v_idx);

    R = min(R + state(:,:,config.b_idx), ones(height, width, 1));
    
    img = cat(3, R, G, B);
    if config.displayLive
        if config.displayZoom ~= 1.0
            img = imresize(img, [height * config.displayZoom, ...
                                 width * config.displayZoom]);
        end
        imshow(img);
    end
    

end

%% representation - recently added by James
% we need a consistent sampling rate: for how long does each data point last?

function strideLength = getStrideLength(phiValue)
    
    % minimum phi value: 0.03
    % maximum phi value: 0.07 (no waves form)

    % if maximum phi value <= 0.06  ->11000
    % if maximum phi value <= 0.055 -> 9000
    % if maximum phi value <= 0.05  -> 8000
    % if maximum phi value <= 0.04  -> 7000
    
    strideLength = 0;
    
    if phiValue <= 0.04
        strideLength = 7000;
    elseif phiValue <= 0.05
        strideLength = 8000;
    elseif phiValue <= 0.055
        strideLength = 9000;
    else 
        strideLength = 11000;
    end
    
end

function [output_data, state]  = runReservoir(state, input_data, maxPhiValue, config)
    
    trace_y = config.trace_y;
    trace_x = config.trace_x;
    u_idx   = config.u_idx;
    
    % convert magnitude of each data point to frequency of waves?
    % let's assume all data is in the range 0..1
    
    % 1 -> maximum frequency (i.e. minimum time between waves)
    % 0 -> minimum frequency (and so will be dictated by our sampling rate)
    
    % need the mapping from phi value to stride length (i.e. how close together can our waves be?
    
    strideLength = getStrideLength(maxPhiValue);
    numSamples   = config.numSamplesPerEpoch;    % at least 3
    timeLength   = strideLength * numSamples;
    
    output_data  = zeros(1, length(input_data) * timeLength); % currently assuming only 1 output
%     output_data  = zeros(config.numOutputs, length(input_data) * timeLength); % currently assuming only 1 output
    
    t = 1;
    for t_i = 1:length(input_data)
        
        mag = input_data(t_i);
        
        if mag > 0
            [quantity, delay] = convertMagToQuant(mag, timeLength, strideLength);
            if quantity > 0
                for q = 1:quantity
                    % stimulate the inputs
                    % wait \strideLength\ timesteps
                    % repeat q times

                    state = stimulateInput(state, 1.0, config);

                    for i = 1:delay
                        state = oreg_step(state, config);
                        
                        if mod(t, config.stride) == 0
                            updateDisplay(state, config);
                        end

                        output_data(t) = state(trace_y, trace_x, u_idx);
                        
                        t = t + 1;
                    end

                end
            else
                for i = 1:timeLength
                    state = oreg_step(state, config);
                    
                    if mod(t, config.stride) == 0
                        updateDisplay(state, config);
                    end
                    
                    output_data(t) = state(trace_y, trace_x, u_idx);
                    t = t + 1;
                end
            end
        else
            for i = 1:timeLength                
                state = oreg_step(state, config);
                if mod(t, config.stride) == 0
                    updateDisplay(state, config);
                end
                
                output_data(t) = state(trace_y, trace_x, u_idx);
                t = t + 1;

            end
        end
        
    end
    
end

function output_data = processOutputs(outputStream, maxPhiValue, config)

    % mag = freq / sampleRate;
    
    % - if we receive a wave every /sampleRate/ timesteps, we have the maximum frequency (i.e. mag=1)
    % - if we receive no waves in /sampleRate/ timesteps, we have the minimum frequency (i.e. mag=0)
    % - counting the number of waves in multiple rounds of /sampleRate/ timesteps, we can produce a linear mapping
    % between frequency of waves and magnitude in the range 0..1

    % convert frequency of waves to magnitude of data
    strideLength = getStrideLength(maxPhiValue); 
    numSamples   = config.numSamplesPerEpoch;    % at least 3
    timeLength   = strideLength * numSamples;
    threshold    = config.waveThreshold;         % how high does the concentration need to be to register as a wave?
    output_data  = zeros(1, ceil(length(outputStream) / timeLength));
    o = 1;
    for t = 1:timeLength:length(outputStream)
        inWave = false;
        waves  = 0;
        for i = 1:timeLength
            % count number of waves
            u = outputStream(t + i - 1);
            if ~inWave && u > threshold
                % rising edge
                inWave = true;
                
                % increment wave counter
                waves = waves + 1;
            end
            
            if inWave && u <= threshold
                % falling edge
                inWave = false;
            end
        end
        output_data(o) = convertQuantToMag(waves, timeLength, strideLength);
        o = o + 1;
    end
    
    
    
end

function mag = convertQuantToMag(quantity, timeLength, strideLength)

    % having measured \quantity\ waves in \timeLength\ timesteps with \strideLength\ stride, what is
    % the underlying magnitude which was inputted?
    
    %  0 -> 0.0
    %  timeLength / strideLength -> 1.0
    %  0..1 -> mag = quant.strideLength / timeLength
    
    mag = quantity * strideLength / timeLength;


end

function [quantity, delay] = convertMagToQuant(mag, timeLength, strideLength)

    % given \mag\ input (0..1), how many waves should we produce in a given length of time?
    % 1.0 and 0.0 are the easiest to calculate: 
    %    0.0 -> 0 waves
    %    1.0 -> timeLength / strideLength waves
    %    for 0..1 : linear mapping between extreme values ^^
    %       -- mag.timeLength / strideLength
    
    freq = mag * timeLength / strideLength; % frequency of waves
    quantity = floor(freq);                 % number of waves to produce
    
    % how many time steps do I need to wait after triggering the wave before I can trigger the next
    % one?
    
    % 1.0 -> strideLength
    % 0.5 -> strideLength * 2
    % 0.0 -> timeLength
    
    % \quantity\ waves to produce in \timeLength\ steps
    % timeLength / freq = time to wait?
    
    % delay = stride / mag (for mag > 0)
    
    delay = timeLength;
    if mag > 0
        delay = strideLength / mag;
        if delay > timeLength
            delay = timeLength;
        end
    end
end


%% step function
function state = oreg_step(state, height, width, config)

    epsilon    = config.epsilon;
    grid_2     = config.grid_2;
    diff_coeff = config.diff_coeff;
    f          = config.f;
    q          = config.q;
    
    u_idx      = config.u_idx;
    v_idx      = config.v_idx;
    b_idx      = config.b_idx;
    p_idx      = config.p_idx;
    
    
    % store changes for t+1 in here
    del_uv     = zeros(height, width, 4);   %need 4 layers here, but we don't change the last 2
    
    for row = 1:height
        cprev    = state(row, 1, u_idx);   % x - 1
        cnext    = state(row, 2, u_idx);   % x + 1
        
        if row == 1, rprev_i = 1; else rprev_i = row - 1; end
%         rprev_i  = iif(row == 1, 1, row - 1);  % boundary condition
        if row == height, rnext_i = height; else rnext_i = row + 1; end
%         rnext_i  = iif(row == height, height, row + 1);  % boundary condition
        
        for col = 1:width
            this_u   = cnext;   % already calculated this value 
            this_v   = state(row, col, v_idx);
            if col == width, cnext_i = width; else cnext_i  = col + 1; end
%             cnext_i  = iif(col == width, width, col + 1);
            cnext    = state(row, cnext_i, u_idx);
            rprev    = state(rprev_i, col, u_idx);
            rnext    = state(rnext_i, col, u_idx);
            phi      = state(row, col, p_idx);       % illumination here
            boundary = state(row, col, b_idx);       % is this a hard boundary / vesicle?
            
            
            if boundary == 1
                del_uv(row, col, u_idx) = -this_u;   % kill the wave
                del_uv(row, col, v_idx) = -this_v;   % 
            else
                laplacian = (rprev + rnext + cprev + cnext - (4 * this_u));
                laplacian = laplacian / grid_2;
                
                del_uv(row, col, u_idx) = ( (...
                                        this_u - (this_u * this_u) - (f * this_v + phi) *  ...
                                        ((this_u - q) / (this_u + q)) ...
                                    ) / epsilon ...
                                    ) + (diff_coeff * laplacian);
                
                del_uv(row, col, v_idx) = this_u - this_v;
                
            end
            
            cprev = this_u;
            
        end
        
    end
    
    state = state + (del_uv * config.speed);
end
