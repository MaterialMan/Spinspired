function[final_states,individual] = collectSTOStates(individual,input_sequence,config,target_output)

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
    
    %if single input entry, add previous state
    if size(input_sequence,1) == 1
        input_sequence = [zeros(size(input_sequence)); input_sequence];
    end
    
    for i= 1:config.num_reservoirs
        if size(input_sequence,1) == 2
            states{i} = individual.last_state{i};
        else
            states{i} = zeros(size(input_sequence,1),individual.nodes(i));
        end
    end

    %states{1} = zeros(size(input_sequence,1),individual.nodes(1) +individual.n_input_units);
    
    % iterate through subreservoirs
    for i = 1:config.num_reservoirs
        
        % change material files
        changeMaterialFile(file_path, individual, i, config);
        
        % change system files
        changeInputFile(file_path, individual, input_sequence, i, config);
        
        % change/update source/input files
        if i > 1
            changeSourceFile(file_path, individual, input_sequence, i, states{i-1},config)
        else
            changeSourceFile(file_path, individual, input_sequence, i, states{i},config)
        end
        

        
        % rewrite files and run
        c1 = strcat('mv "',file_path,'material_file" "',file_path,'material.mat"');
        c2 = strcat('mv "',file_path,'input_file" "',file_path,'input"');
        system(c1); system(c2);
        
        % run vampire!
        command = strcat('cd "', file_path,'" && ./vampire-serial');
        
        [status, cmdout] = system(command);
        
        if status == 0 % if run successful, read reservoir_output file and store in final_states matrix
            output_file = fopen(strcat(file_path,'reservoir_output.txt'), 'r');
            formatSpec = [repmat('%f ',1,individual.nodes(i))];
            size_states = [individual.nodes(i) Inf];
            tmp_states = fscanf(output_file,formatSpec,size_states)';
            fclose(output_file);
            
            tmp_states(isnan(tmp_states)) = 0;
            tmp_states(isinf(tmp_states)) = 0;
            
            % separate 3 directional states
            tmp_x_states = tmp_states(1 : config.interpolation_length*(size(input_sequence,1)) ,:);
            tmp_y_states = tmp_states(config.interpolation_length*(size(input_sequence,1))+1    :...
                config.interpolation_length*(size(input_sequence,1))*2 ,:);
            tmp_z_states = tmp_states(config.interpolation_length*(size(input_sequence,1))*2 +1 : end ,:);
            
            x_states = zeros(size(input_sequence,1),individual.nodes(i));
            y_states = zeros(size(input_sequence,1),individual.nodes(i));
            z_states = zeros(size(input_sequence,1),individual.nodes(i));

            for time = 1:size(input_sequence,1)
                x_states(time,:) = tmp_x_states(config.interpolation_length*time,:);
                y_states(time,:) = tmp_y_states(config.interpolation_length*time,:);
                z_states(time,:) = tmp_z_states(config.interpolation_length*time,:);
%                 z_states(time,:) = 0;
%                 for partial = 0:config.interpolation_length-1
%                     z_states(time,:) = z_states(time,:) + tmp_z_states((config.interpolation_length*time) + partial,:);
%                 end
            end
            
            tmp_states = []; % release
            states{i} = [];
            
            % use only the desired direction(s)
            c.preprocess = 'rescale_diff';
            c.preprocess_shift = [0 1];
            x_states = featureNormailse(x_states,c);
            y_states = featureNormailse(y_states,c);
            z_states = featureNormailse(z_states,c);
            
            for m = 1:length(config.read_mag_direction)
                switch(config.read_mag_direction{m})
                    case 'x'
                        states{i} = [states{i} x_states];
                    case 'y'
                        states{i} = [states{i} y_states];
                    case 'z'
                        states{i} = [states{i} z_states];
                end
            end
            
            % plot states as surface plot
            if config.plot_states
%                 subplot(2,2,2)
%                 plot(input_sequence)
%                 xlim([0 50]) 
%                 title('Input')
                
                subplot(1,2,1)
                plot(tmp_z_states)
                field = individual.applied_field_strength(i);
                title("S_z states, applied field = "+field+"T")
                xlim([25*config.interpolation_length 50*config.interpolation_length]) 
                xlabel('Simulation timesteps')
                ylabel('Z component of the spin from simulation')
                drawnow
                
                subplot(1,2,2)
                plot(z_states)
                title('Scaled S_z states shortened and scaled')
                xlim([25 50]) 
                xlabel('Reservoir output timesteps')
                ylabel('Reservoir output')
                drawnow
                
                
                
%                 if ~config.parallel% surf plot of magnetisation direction
%                     
%                     subplot(2,2,4)
%                     node_grid_size = sqrt(size(x_states,2));
%                     [X,Y] = meshgrid(1:1:node_grid_size, 1:1:node_grid_size);
%                     h = quiver3(X,Y,zeros(size(X)),reshape(x_states(1,1:end),node_grid_size,node_grid_size)...
%                         ,reshape(y_states(1,1:end),node_grid_size,node_grid_size),...
%                         reshape(z_states(1,1:end),node_grid_size,node_grid_size));
%                     
%                     set(gca, 'xDir', 'reverse',...
%                         'camerapositionmode','manual','cameraposition',[1 1 1.25]);
%                     
%                     axis([1 node_grid_size 1 node_grid_size -1.25 1.25]);
%                     zticks([-1 0 1])
%                     zticklabels({min(min(states{i}(:,1:end))) ,0, max(max(states{i}(:,1:end)))})
%                     grid off
%                     
%                     for t = 1:size(states{i},1)
%                         %newH = reshape(states{i}(t,:),node_grid_size,node_grid_size);
%                         set(h,'UData',reshape(x_states(t,:),node_grid_size,node_grid_size));
%                         set(h,'VData',reshape(y_states(t,:),node_grid_size,node_grid_size));
%                         set(h,'WData',reshape(z_states(t,:),node_grid_size,node_grid_size));
%                         drawnow
%                     end
%                 end
            end
        else
            states{i} = zeros(length(input_sequence),config.num_nodes(i));
            fprintf('simulation failed\n')
            cmdout
        end
        
    end
    
    if individual.core_indx ~= 0
    % clean up files
%     rmpath(file_path);
%     
%     try rmdir(file_path,'s');
%           
%     catch
%         warning('Error: Did not remove folder but contiued.\n')
%     end
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


%% change input signal to vampire simulator
function changeSourceFile(file_path, individual, input_sequence, indx, previous_states,config)

input_file = fopen(strcat(file_path,'sourcefield.txt'), 'w');


% write input locations
% inputspos = individual.xy{indx};
% only works on square number of nodes at the moment
inputspos = config.num_nodes;
% 
% inputspos = [inputspos' [-2;-1]]; % print xy coords
% 
% s = [repelem("%3.4f",length(inputspos)-1) '%d' '\n'];
% fprintf(input_file,strjoin(s), inputspos');

% write input values for each location
timesteps = 0:size(input_sequence, 1)-1;

switch(individual.architecture)
    case 'ensemble'
        input = config.input_scaler * (individual.input_scaling(indx)*(individual.input_weights{indx}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')');
        
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
% 
% if config.plot_states
%     subplot(2,2,1)
%     imagesc(rot90(reshape((max(abs(input))),node_grid_size,node_grid_size),3)')    
%     title('Input location')
%     colormap(bluewhitered)
%     colorbar
%     drawnow
% end

% inputfields = [input timesteps'];
inputfields = [input];

interp = config.interpolation_length;
input_length = (size(input_sequence, 1))* interp;
inputfields = zeros(input_length,individual.nodes(indx));

% for time = 0:config.washout
%     for partial = 1:config.interpolation_length
%         inputfields(time*interp+partial,:) = 0.0;
%     end
% end
for time = 0:size(input_sequence, 1)-1
    for partial = 1:config.interpolation_length
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
input_source=fopen(strcat(file_path,'default_input'),'r');
input_file=fopen(strcat(file_path,'input_file'),'w');

while ~feof(input_source)
    
    l=fgetl(input_source); % get line from base file, check if needs to be rewritten

%     % crystal type
%     if contains(l,'create:crystal-structure')
%         switch(individual.material_type{indx})
%             case 'random_alloy'
%                 l = sprintf('create:crystal-structure = fcc'); % should always be fcc crystal
%             otherwise
%                 l = sprintf('create:crystal-structure = %s', config.crystal_structure{indx});
%         end
%     end
%     % periodic boundaries
%     if contains(l,'#create:periodic-boundaries-x') && individual.periodic_boundary(indx,1)
%         l = sprintf('create:periodic-boundaries-x');
%     end
%     if contains(l,'#create:periodic-boundaries-y') && individual.periodic_boundary(indx,2)
%         l = sprintf('create:periodic-boundaries-y');
%     end
%     if contains(l,'#create:periodic-boundaries-z') && individual.periodic_boundary(indx,3)
%         l = sprintf('create:periodic-boundaries-z');
%     end
%     
%     %core shell details
%     if contains(l,'#create:sphere') && config.core_shell(indx)
%         l = sprintf('create:sphere');
%     end
%     if contains(l,'#dimensions:particle-size') && (config.core_shell(indx) || ~strcmp(individual.material_shape{indx},'film'))
%         l = sprintf('dimensions:particle-size = %d !nm',individual.particle_size(indx));
%     end
%     
%     %create shape
%     if contains(l,'#create:sphere') && ~strcmp(individual.material_shape{indx},'film')
%         l = sprintf(strcat('create:',individual.material_shape{indx}));
%     end
    
    % system dimensions
%     if contains(l,'dimensions:unit-cell-size')
%         l = sprintf('dimensions:unit-cell-size = %.3f %s', config.unit_cell_size(indx),config.unit_cell_units{indx});
%     end
    if contains(l,'dimensions:system-size-x')
        l = sprintf('dimensions:system-size-x = %d %s', individual.system_size(indx,1), config.size_units{indx,1});
    end
    if contains(l,'dimensions:system-size-y')
        l = sprintf('dimensions:system-size-y = %d %s', individual.system_size(indx,2), config.size_units{indx,2});
    end
    if contains(l,'dimensions:system-size-z')
        l = sprintf('dimensions:system-size-z = %d %s', individual.system_size(indx,3), config.size_units{indx,3});
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
%         disp(individual.time_steps_increment(indx) * size(input_sequence,1) + individual.time_steps_increment(indx));
        config.total_time_steps = individual.time_steps_increment(indx) * (size(input_sequence,1))* config.interpolation_length;
%         disp(config.total_time_steps);
        l = sprintf('sim:total-time-steps = %d', config.total_time_steps);
    end
    % applied field
    if contains(l,'sim:applied-field-strength') 
        l = sprintf('sim:applied-field-strength = %.2f !T' , individual.applied_field_strength(indx));
    end
     if contains(l,'sim:spin-transfer-torque-polarization-unit-vector') 
        l = sprintf('sim:spin-transfer-torque-polarization-unit-vector = %s' , config.spin_transfer_unit_vector{indx});
    end
    
    % material and unit cell file
    if contains(l,'material:file') 
        l = sprintf('material:file =%s' , config.mat_file);
    end
    if contains(l,'material:unit-cell-file') 
        l = sprintf('material:unit-cell-file = "%s"' , config.unit_cell_file);
    end
    
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
mat_file=fopen(strcat(file_path,'material_file'),'w');

        mat_source=fopen(strcat(file_path,'default_material.txt'),'r');

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
    l = sprintf('material[1]:spin-transfer-relaxation-torque=%s', config.spin_transfer_relaxation_torque{indx});
end

if contains(l,'material[1]:spin-transfer-precession-torque=')
    l = sprintf('material[1]:spin-transfer-precession-torque=%s', config.spin_transfer_precession_torque{indx});
end

if contains(l,'material[1]:spin-transfer-torque-asymmetry=')
    l = sprintf('material[1]:spin-transfer-torque-asymmetry=%s', config.spin_transfer_torque_asymmetry{indx});
end

fprintf(mat_file,'%s \n',l);  % print line to file



fclose(mat_source);
fclose(mat_file);

end

%% change the 
function changeUnitCellFile(file_path, individual, indx, config)


end