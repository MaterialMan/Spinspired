function[final_states,individual] = collectSTOStates(individual,input_sequence,config,target_output)

%% check batch path - for multiple sessions of vampire using multiple cores (stops interference)

%set up the new batch path
if ~isempty(individual.batch_path)
    
    % assign path
    batch_path = individual.batch_path;
    
    %if single input entry, add previous state
    if size(input_sequence,1) == 1
        input_sequence = [zeros(size(input_sequence)); input_sequence];
    end
    
    % add bias node to input
    input_sequence = [input_sequence repmat(config.bias,size(input_sequence,1),1)]; % add bias node to input
    
    % collect states of system - depending on architecture type
    switch(config.architecture)
        case 'ensemble'
            parfor num_res = 1:config.num_res_in_layer(1)
                states{num_res} = getStates(batch_path,individual.layer(1),input_sequence,num_res,num_res,config);
            end
            
        case {'pipeline','pipeline_IA'} % need to finish
            % cycle through layers of pipeline
            for layer = 1:config.num_layers
                if layer > 1 % use initial input for first reservoir
                    previous_states = [zeros(1,size(states{layer-1},2)); states{layer-1}(1:end-1,:)]; %use previous states as input for current layer...
                    
                    if strcmp(config.architecture,'pipeline_IA')
                        %layer_input = [layer_input previous_states];
                        layer_input = [input_sequence previous_states];
                    else
                        layer_input = previous_states;
                    end
                else
                    layer_input = input_sequence;
                end
                
                % collect states of current layer
                tmp_invd = repmat(individual.layer(layer),1,config.num_res_in_layer(layer));
                if config.num_res_in_layer(layer) > 1
                    parfor num_res = 1:config.num_res_in_layer(layer)
                        tmp_invd(num_res).core_indx = num_res;
                        layer_states{num_res} = getStates(batch_path,tmp_invd(num_res),layer_input,num_res,num_res,config);
                    end
                else
                    layer_states{1} = getStates(batch_path,tmp_invd,layer_input,individual.core_indx,1,config);
                end
                
                % concat states for layer to be used for next layer
                states{layer} = [];
                for j = 1:config.num_res_in_layer(layer)
                    states{layer} = [states{layer} layer_states{j}];
                end
                
            end
            
        case 'tree'
            % cycle through layers of pipeline
            for layer = 1:config.num_layers
                if layer > 1 % use initial input for first reservoir
                    previous_states = [zeros(1,size(states{layer-1},2)); states{layer-1}(1:end-1,:)]; %use previous states as input for current layer...
                    
                    if strcmp(config.architecture,'pipeline_IA')
                        %layer_input = [layer_input previous_states];
                        layer_input = [input_sequence previous_states];
                    else
                        layer_input = previous_states;
                    end
                else
                    layer_input = input_sequence;
                end
                
                % collect states of current layer
                tmp_invd = repmat(individual.layer(layer),1,config.num_res_in_layer(layer));
                parfor num_res = 1:config.num_res_in_layer(layer)
                    tmp_invd(num_res).core_indx = num_res
                    layer_states{num_res} = getStates(batch_path,tmp_invd(num_res),layer_input,num_res,num_res,config);
                end
                
                % concat states for layer to be used for next layer
                states{layer} = [];
                for j = 1:config.num_res_in_layer(layer)
                    states{layer} = [states{layer} layer_states{j}];
                end
                
            end
            
        otherwise
            states{1} = getStates(batch_path,individual.layer(1),input_sequence,individual.core_indx,1,config);
    end
    
else
    states{1} = zeros(size(input_sequence,1), config.total_units + individual.n_input_units);
end


% concat all states for output weights
final_states = [];
for layer= 1:config.num_layers
    final_states = [final_states states{layer}];
    
    %assign last state variable
    individual.last_state{layer} = states{layer}(end,:);
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

function states = getStates(batch_path,individual,input_sequence,core_indx,indx,config)

%% find or assign new cores
% check cores on path
current_core = strcat('core',num2str(core_indx));
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
% add to path
%             addpath(genpath(new_dest_path));
%             % copy over files
%             movefile(strcat(new_dest_path,'/',default_core,'.txt'),strcat(new_dest_path,'/',current_core,'.txt'))
%           
    % assign new path
    file_path = new_dest_path;
end

% change material files
changeMaterialFile(file_path, individual, indx, config);

% change system files
changeInputFile(file_path, individual, input_sequence, indx, config);

% update ucf file
changeUnitCellFile(file_path, individual, indx, config)

% change/update source/input files
changeSourceFile(file_path, individual, input_sequence, indx,config)

% rewrite files and run
fclose('all');
c1 = strcat('mv "',file_path,'material_file" "',file_path,'material.mat"');
c2 = strcat('mv "',file_path,'input_file" "',file_path,'input"');
system(c1); system(c2);

% run vampire!
command = strcat('cd "', file_path,'" && ./vampire-STO-serial');

[status, cmdout] = system(command);

if status == 0 % if run successful, read reservoir_output file and store in final_states matrix
    output_file = fopen(strcat(file_path,'reservoir_output.txt'), 'r');
    formatSpec = [repmat('%f ',1,individual.nodes(indx))];
    size_states = [individual.nodes(indx) Inf];
    tmp_states = fscanf(output_file,formatSpec,size_states)';
    fclose(output_file);
    
    if isempty(tmp_states) || (size(tmp_states,1) ~=individual.interpolation_length(indx)*(size(input_sequence,1))*3)
        tmp_states = zeros(individual.interpolation_length(indx)*(size(input_sequence,1))*3,individual.nodes(indx));
    end
    
    tmp_states(isnan(tmp_states)) = 0;
    tmp_states(isinf(tmp_states)) = 0;
    
    % separate 3 directional states
    tmp_x_states = tmp_states(1:individual.interpolation_length(indx)*(size(input_sequence,1)),:);
    tmp_y_states = tmp_states(individual.interpolation_length(indx)*(size(input_sequence,1))+1:...
        individual.interpolation_length(indx)*(size(input_sequence,1))*2,:);
    tmp_z_states = tmp_states(individual.interpolation_length(indx)*(size(input_sequence,1))*2 +1:end,:);
    
    x_states = zeros(size(input_sequence,1),individual.nodes(indx));
    y_states = zeros(size(input_sequence,1),individual.nodes(indx));
    z_states = zeros(size(input_sequence,1),individual.nodes(indx));
    
    % apply envelope
    env_points = 10;
    switch(config.envelope)
        case 'upper'
            [yupper] = envelope(tmp_z_states,env_points,'peak');
            tmp_z_states = [yupper];
        case 'lower'
            [~,ylower] = envelope(tmp_z_states,env_points,'peak');
            tmp_z_states = [ylower];
        case 'both'
            [yupper,ylower] = envelope(tmp_z_states,env_points,'peak');
            tmp_z_states = [yupper,ylower];
        case 'movmean'
            tmp_z_states = movmean(tmp_z_states,[env_points 0],1);
        case 'norm' % Euclidean distance between STO 1 and all others
            tmp_z_states = tmp_z_states-tmp_z_states(:,1);
        otherwise
            
    end
    
    % remove underlying oscillation
    %             Ts = 1;
    %             x = tmp_z_states(201:individual.interpolation_length*config.extra_washout,:);
    %             subplot(1,4,1)
    %             plot(tmp_z_states)
    %
    %             y = fft(x);
    %             fs = 1/Ts;
    %             f = (0:length(y)-1)*fs/length(y);
    %             subplot(1,4,2)
    %             plot(f,abs(y))
    
    %             x = tmp_z_states;
    %             [b,a] = butter(3,0.01);
    %             tmp_z_states = filter(b,a,x);
    
    %             subplot(1,4,3)
    %             plot(dataOut)
    %
    %             subplot(1,4,4)
    %             dataOut = filter(b,a,tmp_z_states);
    %             plot(dataOut)
    
    % interpolate
    for time = 1:size(input_sequence,1)
        switch(config.STO_readout)
            case 'state_average'
                slice = individual.interpolation_length(indx)*(time-1)+1:individual.interpolation_length(indx)*time;
                x_states(time,:) = mean(tmp_x_states(slice,:));
                y_states(time,:) = mean(tmp_y_states(slice,:));
                z_states(time,:) = mean(tmp_z_states(slice,:));
            case 'hilbert_transform'
                slice = individual.interpolation_length(indx)*(time-1)+1:individual.interpolation_length(indx)*time;
                x_states(time,:) = real(hilbert(tmp_x_states(slice,:)));
                y_states(time,:) = real(hilbert(tmp_y_states(slice,:)));
                z_states(time,:) = real(hilbert(tmp_z_states(slice,:)));
            case 'zero_crossing'
                slice = individual.interpolation_length(indx)*(time-1)+1:individual.interpolation_length(indx)*time;
                x_states(time,:) = zeroCrossing(tmp_x_states(slice,:),1);
                y_states(time,:) = zeroCrossing(tmp_y_states(slice,:),1);
                z_states(time,:) = zeroCrossing(tmp_z_states(slice,:),1);
            otherwise
                % last state
                slice = individual.interpolation_length(indx)*time;
                x_states(time,:) = tmp_x_states(slice,:);
                y_states(time,:) = tmp_y_states(slice,:);
                z_states(time,:) = tmp_z_states(slice,:);
        end
        
    end
    
    % use only the desired direction(s)
    states = [];
    for m = 1:length(config.read_mag_direction)
        switch(config.read_mag_direction{m})
            case 'x'
                states = [states x_states];
            case 'y'
                states = [states y_states];
            case 'z'
                states = [states z_states];
        end
    end
    
    % apply envelope
    %[yupper,ylower] = envelope(states{i},10,'peak');
    
    % plot states as surface plot
    if config.plot_states
        set(0,'currentFigure',config.figure_array(1))
        subplot(2,2,1)
        plot(input_sequence)
        title('Input')
        
        subplot(2,2,3)
        plot(tmp_states)
        field = individual.applied_field_strength(indx);
        title("S_z states, applied field = "+field+"T")
        xlabel('Simulation timesteps')
        ylabel('Z component of the spin from simulation')
        drawnow
        
        subplot(2,2,4)
        plot(states)
        title('Scaled S_z states shortened and scaled')
        xlabel('Reservoir output timesteps')
        ylabel('Reservoir output')
        drawnow
        
    end
    tmp_states = []; % release
else
    states = zeros(size(input_sequence,1),config.num_nodes(indx));
    fprintf('simulation failed\n')
    cmdout
end

% get leak states
if config.leak_on
    states = getLeakStates(states,individual,1);
end

if core_indx ~= 0
    % clean up files
    rmpath(file_path);
    
    try rmdir(file_path,'s');
        
    catch
        warning('Error: Did not remove folder but contiued.\n')
    end
end

end


%% change input signal to vampire simulator
function changeSourceFile(file_path, individual, input_sequence, indx,config)

input_file = fopen(strcat(file_path,'sourcefield.txt'), 'w');


% write input locations
% inputspos = individual.xy{indx};
% only works on square number of nodes at the moment
inputspos = individual.nodes(indx);
%input = config.input_scaler * (individual.input_scaling(indx)*(individual.input_weights{indx}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')');
input = config.input_scaler * (individual.input_scaling(indx)*(input_sequence*individual.input_weights{indx}));        
%
if config.plot_states
    set(0,'currentFigure',config.figure_array(1))
    node_grid_size = sqrt(individual.nodes(indx));
    subplot(2,2,2)
    imagesc(rot90(reshape((max(abs(input))),node_grid_size,node_grid_size),3)')
    title('Input location')
    colormap(bluewhitered)
    colorbar
    drawnow
end

% inputfields = [input timesteps'];
inputfields = [input];

interp = individual.interpolation_length(indx);
input_length = (size(input_sequence, 1))* interp;
inputfields = zeros(input_length,individual.nodes(indx));

for time = 0:size(input_sequence, 1)-1
    for partial = 1:individual.interpolation_length
        inputfields(time*interp+partial,:) = input(time+1,:);
    end
end

a = [input_length, individual.nodes(indx)];
fprintf(input_file,'%d %d \n',a');

s = [repelem("%3.4f",inputspos) '\n'];

fprintf(input_file,strjoin(s), inputfields');

fclose(input_file);
end

%% change input file to vampire simulation
function changeInputFile(file_path, individual, input_sequence, indx, config)

% change parameters in input file
input_source=fopen(strcat(file_path,'default_STO_input'),'r');
input_file=fopen(strcat(file_path,'input_file'),'w');

while ~feof(input_source)
    
    l=fgetl(input_source); % get line from base file, check if needs to be rewritten
    
    % system dimensions
    %     if contains(l,'dimensions:unit-cell-size')
    %         l = sprintf('dimensions:unit-cell-size = %.3f %s', config.unit_cell_size(indx),config.unit_cell_units{indx});
    %     end
    if contains(l,'dimensions:system-size-x')
        l = sprintf('dimensions:system-size-x = %d %s', individual.system_size(indx,1), config.size_units{1});
    end
    if contains(l,'dimensions:system-size-y')
        l = sprintf('dimensions:system-size-y = %d %s', individual.system_size(indx,2), config.size_units{2});
    end
    if contains(l,'dimensions:system-size-z')
        l = sprintf('dimensions:system-size-z = %d %s', individual.system_size(indx,3), config.size_units{3});
    end
    
    % sim details
    if contains(l,'sim:temperature')
        l = sprintf('sim:temperature = %d', individual.temperature(indx));
    end
    if contains(l,'sim:time-step')
        if contains(l,'sim:time-steps-increment')
            l = sprintf('sim:time-steps-increment = %d', individual.time_steps_increment(indx));
        else
            l = sprintf('sim:time-step = %d %s', config.time_step,config.time_units);
        end
    end
    if contains(l,'sim:total-time-steps')
        config.total_time_steps = individual.time_steps_increment(indx) * (size(input_sequence,1))* individual.interpolation_length(indx);
        l = sprintf('sim:total-time-steps = %d', config.total_time_steps);
    end
    % applied field
    if contains(l,'sim:applied-field-strength')
        l = sprintf('sim:applied-field-strength = %.2f !T' , individual.applied_field_strength(indx));
    end
    if contains(l,'sim:spin-transfer-torque-polarization-unit-vector')
        l = sprintf('sim:spin-transfer-torque-polarization-unit-vector = %s' , config.spin_transfer_unit_vector{1});
    end
    
    % material and unit cell file
    if contains(l,'material:file')
        l = sprintf('material:file =%s' , config.mat_file);
    end
    %     if contains(l,'material:unit-cell-file')
    %         l = sprintf('material:unit-cell-file = "%s"' , config.unit_cell_file);
    %     end
    
    % plot system
    %     if contains(l,'sim:cells-source-output') % plot magnetic moments
    %         if config.plt_system
    %             l = sprintf('sim:cells-source-output = true');
    %         else
    %             l = sprintf('#sim:cells-source-output = true');
    %         end
    %     end
    %     if contains(l,'config:macro-cells') % plot magnetic moments
    %         if contains(l,'config:macro-cells-output-rate') % plot magnetic moments
    %             if config.plt_system
    %                 l = sprintf('config:macro-cells-output-rate = %d', config.plot_rate);
    %             else
    %                 l = sprintf('#config:macro-cells-output-rate = 1');
    %             end
    %         else
    %             if config.plt_system
    %                 l = sprintf('config:macro-cells');
    %             else
    %                 l = sprintf('#config:macro-cells');
    %             end
    %         end
    %     end
    
    fprintf(input_file,'%s \n',l);  % print line to file
end

fclose(input_source);
fclose(input_file);
end

%% change material file to vampire simulation
function changeMaterialFile(file_path, individual, indx, config)

% change parameters in material file
mat_source=fopen(strcat(file_path,'default_STO_material.txt'),'r');
mat_file=fopen(strcat(file_path,'material_file'),'w');

while ~feof(mat_source)
    
    l=fgetl(mat_source); % get line from base file, check if needs to be rewritten
    
    if contains(l,'material[1]:damping-constant=')
        l = sprintf('material[1]:damping-constant=%.3f', individual.damping(indx));
    end
    if contains(l,'material[1]:exchange-matrix[1]')
        l = sprintf('material[1]:exchange-matrix[1]=%s', individual.exchange(indx));
    end
    if contains(l,'material[1]:atomic-spin-moment')
        l = sprintf('material[1]:atomic-spin-moment=%.3f !muB', individual.magmoment(indx));
    end
    if contains(l,'material[1]:second-order-uniaxial-anisotropy-constant=')
        l = sprintf('material[1]:second-order-uniaxial-anisotropy-constant=%s', individual.anisotropy(indx));
    end
    
    if contains(l,'material[1]:initial-spin-direction=')
        l = sprintf('material[1]:initial-spin-direction=%s', config.initial_spin_direction{indx});
    end
    
    
    % some STO things that were not there for MM
    if contains(l,'material[1]:spin-transfer-relaxation-torque=')
        l = sprintf('material[1]:spin-transfer-relaxation-torque=%s', individual.spin_transfer_relaxation_torque(indx));
    end
    
    if contains(l,'material[1]:spin-transfer-precession-torque=')
        l = sprintf('material[1]:spin-transfer-precession-torque=%s', individual.spin_transfer_precession_torque(indx));
    end
    
    if contains(l,'material[1]:spin-transfer-torque-asymmetry=')
        l = sprintf('material[1]:spin-transfer-torque-asymmetry=%s', individual.spin_transfer_torque_asymmetry(indx));
    end
    
    fprintf(mat_file,'%s \n',l);  % print line to file
end

fclose(mat_source);
fclose(mat_file);

end

%% change the
function changeUnitCellFile(file_path, individual, indx, config)
% change parameters in material file

ucf_file=fopen(strcat(file_path,'STO-array.ucf'),'w');

% set unit cell size
fprintf(ucf_file,'# unit cell size \n');
fprintf(ucf_file,'%d %d %d \n',individual.unit_cell_size(indx,:));

% unit cell vectors
fprintf(ucf_file,'# unit cell vectors \n');
fprintf(ucf_file,'%d %d %d \n',eye(3));

% set atom array
[x,y] = meshgrid(0:1/sqrt(individual.nodes(indx)):(1-1/sqrt(individual.nodes(indx))));
atoms = [((1:individual.nodes(indx))-1)' x(:) y(:)];
fprintf(ucf_file,'# Atoms \n');
fprintf(ucf_file,strcat(num2str(individual.nodes(indx))," ",num2str(1),'\n'));

fprintf(ucf_file,'%d %.3f %.3f %d %d %d %d\n',[atoms'; zeros(2,size(atoms,1)); ones(2,size(atoms,1))]);

% define interactions
G = digraph(full(individual.W{1}));
fprintf(ucf_file,'# interactions \n');
fprintf(ucf_file, strcat(num2str(size(individual.xy{indx},1))," ",'normalised-isotropic \n'));
%fprintf(ucf_file,'%d %d %d %d %d %.3f %.3f\n',[(1:size(individual.xy{indx},1))' individual.xy{indx} zeros(size(individual.xy{indx},1),3) G.Edges.Weight.*individual.W_scaling(indx,indx)]');

[p,i] = find(individual.W{1});
xy_weights = zeros(length(individual.xy{indx}),1);
ind_w = [p i]-1;
for c = 1:size(ind_w,1)
    %xy_weights(sum(ismember(individual.xy{indx},g(c,:)),2) == 2) = individual.W{1}(g(c,1)+1,g(c,2)+1);
    xy_weights(ismember(individual.xy{indx},ind_w(c,:),'rows')) = individual.W{1}(ind_w(c,1)+1,ind_w(c,2)+1);
end
tmp_print = [(1:size(individual.xy{indx},1))' individual.xy{indx} zeros(size(individual.xy{indx},1),3) xy_weights.*individual.W_scaling(indx)]';
fprintf(ucf_file,'%d %d %d %d %d %d %.3f\n',tmp_print);

% close file
fclose(ucf_file);

end