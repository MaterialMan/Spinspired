clear
close all

%set random seed for experiments
rng(1,'twister');

%% Setup
config.parallel = 0;                        % use parallel toolbox

%start paralllel pool if empty
if isempty(gcp) && config.parallel %
    parpool('local',4,'IdleTimeout', Inf); % create parallel pool
end

% type of network to evolve
config.res_type = 'multiMM';                % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = {[196]};                  % num of nodes in each sub-reservoir, e.g. if config.num_nodes = {10,5,15}, there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. For one reservoir, sate as a non-cell, e.g. config.num_nodes = 25
config = selectReservoirType(config);   % collect function pointers for the selected reservoir type

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = 'test_pulse';          % Task to evolve for
config.figure_array = [figure figure figure];
config.preprocess ='';

config.metrics = {'linearMC'};

%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = 1;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

%% sweep film sizes
config.test = 1;

%% damping effect
damping = 0.1;%linspace(0.1,1,4);

% manual quad list
quad_list = [1 3 34 36;...
    7 9 92 94;...
    9 13 184 188];

% define list of experiments
input_list = {'c','lu','ld','ru','rd'};

input_idx = 1;
width  = 14;
height = 14;

config.num_nodes{1} = width.^2;

% get any additional params stored in getDataSetInfo.m. This might include:
% details on reservoir structure, extra task variables, etc.
[config] = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

config.add_input_states = 0;

default_pop = config.createFcn(config);

% set dimensions
tp_dim = height/width;
dim_list = [1;tp_dim];

% make rectangle
coord = createRect(dim_list(2),dim_list(1));
x = coord(:,1);
y = coord(:,2);

% % reset input weights
square_film_dimensions = width;
[xq,yq] = meshgrid(linspace(0,1,square_film_dimensions),linspace(0,1,square_film_dimensions));
[in,on] = inpolygon(xq,yq,x,y);
inputs_in_use = find(in | on);

% define input
% grad quadrants
loc_mat = reshape(inputs_in_use,height,width);

switch input_list{input_idx}
    case 'c'
        if mod(square_film_dimensions,2) == 0
            %left (up & down)
            lu = round(length(inputs_in_use)/2)-floor(height/2);
            % make sure input is symmetrical in film
            if mod(height,2) == 0
                input_loc = inputs_in_use([lu lu+1 lu+height lu+height+1]); %[lu ld ru rd]
            else
                input_loc = inputs_in_use([lu lu+height]);
            end
        else
            lu = round(length(inputs_in_use)/2);
            % make sure input is symmetrical in film
            input_loc = inputs_in_use([lu lu+1]); %[lu ld ru rd]
        end
        
    case {'lu','ld','ru','rd'}
        quad_height = size(loc_mat,1);
        quad_width = size(loc_mat,2);
        
        try C = mat2cell(loc_mat,[floor(quad_height/2) floor(quad_height/2)],[ceil(quad_width/2) ceil(quad_width/2)]);
            
            if mod(size(C{input_idx-1},2),2) == 0
                input_loc = C{input_idx-1}(round(numel(C{input_idx-1})/2) - floor(size(C{input_idx-1},1)/2));
            else
                input_loc = C{input_idx-1}(round(numel(C{input_idx-1})/2));
            end
            
        catch
            switch length(inputs_in_use)
                case 36
                    input_loc = inputs_in_use(quad_list(1,input_idx));
                case 100
                    input_loc = inputs_in_use(quad_list(2,input_idx));
                case 196
                    input_loc = inputs_in_use(quad_list(3,input_idx));
            end
        end
end

% stop broadcast variables
geo_width = dim_list(2);
geo_height = dim_list(1);

figure(config.figure_array(1))
s = zeros(square_film_dimensions);
s(inputs_in_use) = -1;
s(input_loc) = 1;
imagesc(s)
drawnow

%ppm = ParforProgMon('Initial population: ', length(damping));
parfor pop_indx = 1:length(damping)
    warning('off','all')
    
    population(pop_indx) = default_pop;
    population(pop_indx).layer(1).input_scaling = 1;
    population(pop_indx).layer(1).leak_rate = 1;
    population(pop_indx).layer(1).input_weights{1} = zeros(2,config.num_nodes{1});
    
    population(pop_indx).layer(1).input_weights{1}(1,input_loc) = 1;
    
    population(pop_indx).layer(1).geo_width = geo_width;
    population(pop_indx).layer(1).geo_height = geo_height;
    
    population(pop_indx).layer(1).core_indx = pop_indx;
    population(pop_indx).layer(1).damping = damping(pop_indx);
    
    states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
end

%% input pulse
%[X,Y] = meshgrid(1:sqrt(config.num_nodes{1}),1:sqrt(config.num_nodes{1}));
[X,Y] = meshgrid(1:width,1:height);
figure
set(gcf,'color','w')
d = 1;
list = 11:3:20;
max_scaler = max(max(states{d}));
min_scaler = min(min(states{d}));
tl = tiledlayout(2,4);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
for t = 1:4
    %subplot(2,4,t)
    h(t) = nexttile(t);
    surf(X,Y,flip(reshape(states{d}(list(t),inputs_in_use),height,width))...
        ,'FaceAlpha',0.5,'FaceColor','interp');
    set(gca,'visible','off')
    %view([2 2 2])
    zlim([min_scaler max_scaler])
    caxis([min_scaler max_scaler])
    if t ==2
        %sgtitle('Side view')
    end
    set(findall(gca,'type','text'),'visible','on')
    set(gca, 'FontSize',16,'FontName','Arial')
    
    %subplot(2,4,t+4)
    h(t+4) = nexttile(t+4);
    surf(X,Y,flip(reshape(states{d}(list(t),inputs_in_use),height,width))...
        ,'FaceAlpha',0.5,'FaceColor','interp');
    view(2)
    caxis([min_scaler max_scaler])
    set(gca,'visible','off')
    if t ==2
        %sgtitle('Top down view')
    end
    xlabel(strcat('n = ',num2str(list(t))))
    set(findall(gca,'type','text'),'visible','on')
    set(gca, 'FontSize',16,'FontName','Arial')
    %colorbar
    colormap(jet)
end
title(tl,'')
colorbar(tl)

if strcmp(config.dataset,'test_pulse')
    print(strcat('input_response_',num2str(width),'_by_',num2str(height),'_pulse_',input_list{input_idx}),'-dpdf','-bestfit')   
else
    print(strcat('input_response_',num2str(width),'_by_',num2str(height),'_',num2str(config.freq),'Hz_',input_list{input_idx}),'-dpdf','-bestfit')
end

figure
for i = 1:length(inputs_in_use)
    %if sum(states{1}(:,i)) > 0
        Y = fft(states{1}(:,inputs_in_use(i))); %*1e-20
        L = size(states{1}(:,inputs_in_use(i)),1);
        Fs = 1e3;
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        plot(f,P1,'b')
        hold on
    %end
end
hold off
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f(GHz)')
ylabel('|P1(f)|')
ylim([0 0.4])
input_list2 = {'c','Q1','Q2','Q3','Q4'};
title(strcat('Periodogram: cells=',num2str(width),'x',num2str(height),', input=',input_list2{input_idx}))
set(gca, 'FontSize',16,'FontName','Arial')

figure(config.figure_array(3))
z = sum(abs(states{1}(:,inputs_in_use)));
z = reshape(z,height,width);
%z = reshape(z,width,height);
subplot(1,2,1)
surf(X,Y,flip(z)...
    ,'FaceAlpha',0.5,'FaceColor','interp');
view(-10.5,67.3)
colormap(jet)
set(gca,'visible','off')

subplot(1,2,2)
surf(X,Y,flip(z)...
    ,'FaceAlpha',0.5,'FaceColor','interp');
view(2)
%imagesc(z)
colormap(jet)
set(gca,'visible','off')
set(gcf,'PaperOrientation','portrait');
print(strcat('jet_3Dexample_sine_',num2str(width),'_by_',num2str(height),'_',num2str(config.freq),'Hz_',input_list{input_idx}),'-dpdf','-bestfit')

figure(config.figure_array(1))
for i = 1:length(damping)
   subplot(2,2,i)
   plot(states{i});
end



% discrete fourier transform(dft)
figure
for i = 1:length(inputs_in_use)
    %if sum(states{1}(:,i)) > 0
        
        subplot(1,2,1)
       % n = 512;
        y = fft(states{1}(:,inputs_in_use(i)));
        m = abs(y);
        f = (0:length(y)-1)*Fs/length(y); 
        hold on
        plot(f,m)
        
        
        subplot(1,2,2)
        y(m<1e-6) = 0;
        p = unwrap(angle(y)); 
        hold on
        plot(f,p*180/pi)
    %end
end
hold off
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f(GHz)')
ylabel('DFT amplitude')

figure
% Fs = 1000;
% t = 0:1/Fs:1-1/Fs;
x =reshape(states{1}(:,inputs_in_use),size(states{1},1),length(inputs_in_use)); %states{1}; width*height
% [Pxx,F] = periodogram(x,[],length(x),Fs);
% plot(F,10*log10(Pxx))
periodogram(x)
input_list2 = {'c','Q1','Q2','Q3','Q4'};
title(strcat('Periodogram: cells=',num2str(width),'x',num2str(height),', input=',input_list2{input_idx}))
set(gca, 'FontSize',18,'FontName','Arial')
ylim([-120 20])


print(strcat('periodogram_pulse_',num2str(width),'_by_',num2str(height),'_',num2str(config.freq),'Hz_',input_list{input_idx}),'-dpdf','-bestfit')