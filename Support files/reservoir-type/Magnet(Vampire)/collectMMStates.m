function[final_states,individual] = collectMMStates(individual,input_sequence,config,target_output)

%% check batch path - for multiple sessions of vampire using multiple cores (stops interference)
% add batch to path
current_batch = strcat('batch',num2str(individual.batch_num),'.txt');
batch_path = which(current_batch);
batch_path = batch_path(1:end-length(current_batch)-1);

if isempty(batch_path) %&&  individual.core_indx == 1 % create new batch directory
    default_batch = 'batch0.txt'; % default batch
    default_batch_path = which(default_batch);
    path_to_copy = default_batch_path(1:end-length(default_batch)-1);
    new_dest_path = strcat(path_to_copy(1:end-length('batch0')),current_batch(1:end-4));
    
    % make new directory
    copyfile(path_to_copy,new_dest_path);
    % copy over files
    movefile(strcat(new_dest_path,'/',default_batch),strcat(new_dest_path,'/',current_batch))
    % add to path
    addpath(genpath(batch_path));
    
    % assign new path
    batch_path = new_dest_path;    
end

if ~isempty(batch_path)
    
    %% find or assign new cores
    % check cores on path
    current_core = strcat('core',num2str(individual.core_indx));
    core_path = strcat(batch_path,'/Cores/',current_core,'/',current_core,'.txt');
    core_check = exist(core_path,'file'); % will be 7 if folder exists
    
    if core_check == 7
        % exists
        file_path = core_path(1:end-(length(current_core)+4));
    else
        % does not exist
        default_core = 'core0'; % default batch
        path_to_copy = strcat(batch_path,'/Cores/',default_core);
        new_dest_path = core_path(1:end-(length(current_core)+4));
        
        % make new directory
        copyfile(path_to_copy,new_dest_path);
        % copy over files
        movefile(strcat(new_dest_path,'/',default_core,'.txt'),strcat(new_dest_path,'/',current_core,'.txt'))
        % add to path
        addpath(genpath(new_dest_path));
        % assign new path
        file_path = new_dest_path;
    end
    
    states{1} = zeros(size(input_sequence,1),individual.nodes(1) +individual.n_input_units);
    
    % iterate through subreservoirs
    for i = 1:config.num_reservoirs
        
        % change/update files
        if i > 1
            changeSourceFile(file_path, individual, input_sequence, i, states{i-1},config)
        else
            changeSourceFile(file_path, individual, input_sequence, i, states{i},config)
        end
        changeCoFile(file_path, individual, i);
        
        changeInputFile(file_path, individual, input_sequence, config, i);
        
        % rewrite files and run
        c1 = strcat('mv "',file_path,'Co_temp" "',file_path,'Co.mat"');
        c2 = strcat('mv "',file_path,'inputtemp" "',file_path,'input"');
        system(c1); system(c2);
        
        % run vampire!
        command = strcat('cd "', file_path,'" && ./vampire-serial');
        
        [status, ~] = system(command);
        
        if status == 0 % if run successful, read reservoir_output file and store in final_states matrix
            output_file = fopen(strcat(file_path,'reservoir_output.txt'), 'r');
            formatSpec = [repmat('%f ',1,individual.nodes(i))];
            size_states = [individual.nodes(i) Inf];
            tmp_states = fscanf(output_file,formatSpec,size_states)';
            fclose(output_file);
            
            tmp_states(isnan(tmp_states)) = 0;
            tmp_states(isinf(tmp_states)) = 0;
            
            % separate 3 directional states
            x_states = tmp_states(1:size(input_sequence,1),:);
            y_states = tmp_states(size(input_sequence,1)+1:size(input_sequence,1)*2,:);
            z_states = tmp_states(size(input_sequence,1)*2 +1:end,:);
            
            tmp_states = []; % release
            states{i} = [];
            
            % use only the desired direction(s)
            x_states = featureNormailse(x_states,'rescale');
            y_states = featureNormailse(y_states,'rescale');
            z_states = featureNormailse(z_states,'rescale');
            
            for m = 1:length(config.read_mag_direction)
                switch(config.read_mag_direction{m})
                    case 'x'
                        % x_states = featureNormailse(x_states,'rescale');
                        states{i} = [states{i} x_states];
                        %x_states = []; % release
                    case 'y'
                        %y_states = featureNormailse(y_states,'rescale');
                        states{i} = [states{i} y_states];
                        %y_states = []; % release
                    case 'z'
                        %z_states = featureNormailse(z_states,'rescale');
                        states{i} = [states{i} z_states];
                        %z_states = []; % release
                end
            end
            
            % plot states as surface plot
            if config.plot_states
                subplot(2,2,2)
                plot(input_sequence)
                title('Input')
                
                subplot(2,2,3)
                plot(states{i})
                title('Scaled states')
                drawnow
                
                if ~config.parallel% surf plot of magnetisation direction
                    
                    subplot(2,2,4)
                    node_grid_size = sqrt(size(x_states,2));
                    [X,Y] = meshgrid(1:1:node_grid_size, 1:1:node_grid_size);
                    h = quiver3(X,Y,zeros(10),reshape(x_states(1,1:end),node_grid_size,node_grid_size)...
                        ,reshape(y_states(1,1:end),node_grid_size,node_grid_size),...
                        reshape(z_states(1,1:end),node_grid_size,node_grid_size));
                    
                    set(gca, 'xDir', 'reverse',...
                        'camerapositionmode','manual','cameraposition',[1 1 1.25]);
                    
                    axis([1 node_grid_size 1 node_grid_size -1.25 1.25]);
                    zticks([-1 0 1])
                    zticklabels({min(min(states{i}(:,1:end))) ,0, max(max(states{i}(:,1:end)))})
                    grid off
                    
                    for t = 1:size(states{i},1)
                        %newH = reshape(states{i}(t,:),node_grid_size,node_grid_size);
                        set(h,'UData',reshape(x_states(t,:),node_grid_size,node_grid_size));
                        set(h,'VData',reshape(y_states(t,:),node_grid_size,node_grid_size));
                        set(h,'WData',reshape(z_states(t,:),node_grid_size,node_grid_size));
                        drawnow
                    end
                end
            end
        else
            states{i} = zeros(length(input_sequence),config.num_nodes(i));
        end
        
    end
    
    if individual.core_indx ~= 0
    % clean up files
    rmpath(file_path);
    
    try rmdir(file_path,'s');
        
    catch
        warning('Error: Did not remove folder but contiued.\n')
    end
    end

else
    states{1} = zeros(size(input_sequence,1),individual.nodes(1) +individual.n_input_units); 
end

% get leak states
if config.leak_on
    states = getLeakStates(states,individual,input_sequence,config);
end

% concat all states for output weights
final_states = [];
for i= 1:config.num_reservoirs
    final_states = [final_states states{i}];
    
    %assign last state variable
    individual.last_state{i} = states{i}(end,:);
end

% concat input states
if config.add_input_states == 1
    final_states = [final_states input_sequence];
end

if size(input_sequence,1) == 2
    final_states = final_states(end,:); % remove washout
else
    final_states = final_states(config.wash_out+1:end,:); % remove washout
end


end



function changeSourceFile(file_path, individual, input_sequence, indx, previous_states,config)

input_file = fopen(strcat(file_path,'sourcefield.txt'), 'w');

% write input locations
inputspos = individual.xy{indx};

inputspos = [inputspos' [-2;-1]]; % print xy coords

s = [repelem("%3.4f",length(inputspos)-1) '%d' '\n'];
fprintf(input_file,strjoin(s), inputspos');

% write input values for each location
timesteps = 0:size(input_sequence, 1)-1;

switch(individual.architecture)
    case 'ensemble'
        input = individual.input_scaling(indx)*(individual.input_weights{indx}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')';
        
    case 'pipeline_IA'
        if indx < 2
            input = individual.input_scaling(indx)*(individual.input_weights{indx}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')';
        else
            input = previous_states(:,1:end-individual.n_input_units*(size(previous_states,2) > individual.nodes(indx))) * (individual.input_scaling(indx)*individual.input_weights{indx})...
                + individual.W{indx-1,indx}*individual.W_scaling(indx-1,indx)*previous_states;
        end
    case 'pipeline'
        if indx < 2
            input = individual.input_scaling(indx)*(individual.input_weights{indx}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')';
        else
            input =  previous_states(:,1:end-individual.n_input_units*(size(previous_states,2) > individual.nodes(indx))) * (individual.W_scaling(indx-1,indx)*individual.W{indx-1,indx});
        end
end

% change input widths
node_grid_size = sqrt(individual.nodes(indx));
for n = 1:size(input,1)
    m = reshape(input(n,:),node_grid_size,node_grid_size);
    f_pos = find(m);
    input_matrix_2d = m;
    for p = 1:length(f_pos)
        t = zeros(size(m));
        t(f_pos(p)) = m(f_pos(p));
        [t] = adjustInputShape(t,individual.input_widths{indx}(f_pos(p)));
        input_matrix_2d = input_matrix_2d + t;
    end
    input(n,:) = input_matrix_2d(:);
end

if config.plot_states
    subplot(2,2,1)
    imagesc(reshape(max(abs(input)),node_grid_size,node_grid_size))
    title('Input location')
    colormap(bluewhitered)
    colorbar
    drawnow
end

inputfields = [input timesteps'];

s = [repelem("%3.4f",length(inputspos)-1) '%d' '\n'];
fprintf(input_file,strjoin(s), inputfields');

fclose(input_file);
end

function changeInputFile(file_path, individual, input_sequence, config, indx)

% change parameters in input file
base_input=fopen(strcat(file_path,'base_input'),'r');
inputtemp=fopen(strcat(file_path,'inputtemp'),'w');

while ~feof(base_input)
    l=fgetl(base_input); % get line from base file, check if needs to be rewritten
    
    % system parameters
    if contains(l,'dimensions:unit-cell-size')
        l = sprintf('dimensions:unit-cell-size = %.2f !A', config.unitcell_size);
    end
    if contains(l,'dimensions:system-size-x')
        l = sprintf('dimensions:system-size-x = %d %s', individual.system_size(indx), config.size_units);
    end
    if contains(l,'dimensions:system-size-y')
        l = sprintf('dimensions:system-size-y = %d %s', individual.system_size(indx), config.size_units);
    end
    %     if contains(l,'dimensions:system-size-z')
    %         l = sprintf('dimensions:system-size-z = %d %s', individual.system_size(indx), config.size_units);
    %     end
    if contains(l,'cells:macro-cell-size')
        if ~contains(l,'#cells:macro-cell-size')
            l = sprintf('cells:macro-cell-size = %d %s', config.macro_cell_size, config.size_units);
        end
    end
    if contains(l,'create:crystal-structure')
        l = sprintf('create:crystal-structure = %s', config.crystal_structure);
    end
    
    
    % sim details
    if contains(l,'sim:temperature')
        l = sprintf('sim:temperature = %d', individual.temperature(indx));
    end
    if contains(l,'sim:time-step')
        if contains(l,'sim:time-steps-increment')
            l = sprintf('sim:time-steps-increment = %d', config.time_steps_increment);
        else
            l = sprintf('sim:time-step = %d %s', config.time_step,config.time_units);
        end
    end
    if contains(l,'sim:total-time-steps') % MD added
        config.total_time_steps = config.time_steps_increment * length(input_sequence)+ config.time_steps_increment;
        l = sprintf('sim:total-time-steps = %d', config.total_time_steps);
    end
    if contains(l,'sim:applied-field-strength') % MD added
        l = sprintf('sim:applied-field-strength = %.2f !T' , individual.applied_field_strength(indx));
    end
    
    
    % plot system
    if contains(l,'sim:cells-source-output') % plot magnetic moments
        if config.plt_system
            l = sprintf('sim:cells-source-output = true');
        else
            l = sprintf('#sim:cells-source-output = true');
        end
    end
    if contains(l,'config:macro-cells') % plot magnetic moments
        if contains(l,'config:macro-cells-output-rate') % plot magnetic moments
            if config.plt_system
                l = sprintf('config:macro-cells-output-rate = %d', config.plot_rate);
            else
                l = sprintf('#config:macro-cells-output-rate = 1');
            end
        else
            if config.plt_system
                l = sprintf('config:macro-cells');
            else
                l = sprintf('#config:macro-cells');
            end
        end
    end
    
    fprintf(inputtemp,'%s \n',l);  % print line to file
end

fclose(base_input);
fclose(inputtemp);
end

function changeCoFile(file_path, individual, indx)

% change parameters in Co.mat
comat=fopen(strcat(file_path,'Co.txt'),'r');
cotemp=fopen(strcat(file_path,'Co_temp'),'w');

while ~feof(comat)
    l=fgetl(comat); % get line from base file, check if needs to be rewritten
    if contains(l,'material[1]:damping-constant=')
        l = sprintf('material[1]:damping-constant=%d', individual.damping(indx));
    end
    if contains(l,'material[1]:second-order-uniaxial-anisotropy-constant=')
        l = sprintf('material[1]:second-order-uniaxial-anisotropy-constant=%s', individual.anisotropy(indx));
    end
    if contains(l,'material[1]:exchange-matrix[1]')
        l = sprintf('material[1]:exchange-matrix[1]=%s', individual.exchange(indx));
    end
    
    if contains(l,'material[1]:atomic-spin-moment')
        l = sprintf('material[1]:atomic-spin-moment=%d !muB', individual.magmoment(indx));
    end
    if contains(l,'material[1]:material-element') % MD added
        l = sprintf('material[1]:material-element=%s', individual.material_element{indx});
    end
    
    fprintf(cotemp,'%s \n',l);  % print line to file
end

fclose(comat);
fclose(cotemp);
end