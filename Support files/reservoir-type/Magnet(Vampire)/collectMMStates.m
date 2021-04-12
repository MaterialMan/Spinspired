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
    
    % collect states of system - depending on architecture type
    switch(config.architecture)
        case 'ensemble'
            parfor i = 1:config.num_reservoirs % cyclce through all subres - all independent from eachother
                states{i} = getStates(batch_path,individual,input_sequence,states{i},i,config)
            end
        case 'pipeline' % need to finish
            % cycle through layers of pipeline
            for i = 1:config.num_layers
                if i == 1 % use iniital input for first reservoir
                    previous_states = states{i};
                else
                    previous_states = states{i-1}; %use previous states as input for current layer... 
                    % need to add weighting mechanism
                end
                % collect states of current layer
                parfor j = 1:individual.num_reservoirs_in_layer(i)
                        layer_states{j} = getStates(batch_path,individual(i),input_sequence,previous_states,j,config);
                end
                % concat states for layer to be used for next layer
                for j = 1:config.num_reservoirs
                    states{i} = [states{i} layer_states{j}];
                end
            end
            
        otherwise
            states{1} = getStates(batch_path,individual,input_sequence,states{1},1,config);
    end

else
    states{1} = zeros(size(input_sequence,1),sum(config.num_nodes) +individual.n_input_units);
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


%% assess ensemble reservoirs
function states = getStates(batch_path,individual,input_sequence,previous_states,i,config)
    % iterate through subreservoirs
   % parfor i = 1:config.num_reservoirs
        
        %% find or assign new cores
        % check cores on path
        current_core = strcat('core',num2str(individual.core_indx(i)));
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
        
        % change/update source/input files
        changeSourceFile(file_path, individual, input_sequence, i, previous_states,config)
        
        % change material files
        changeMaterialFile(file_path, individual, i, config);
        
        % change system files
        changeInputFile(file_path, individual, input_sequence, i, config);
        
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
            
            tmp_states = tmp_states.*1e20;
            
            % separate 3 directional states
            x_states = tmp_states(1:size(input_sequence,1),:);
            y_states = tmp_states(size(input_sequence,1)+1:size(input_sequence,1)*2,:);
            z_states = tmp_states(size(input_sequence,1)*2 +1:end,:);
            
            tmp_states = []; % release
            states = [];
            
            % use only the desired direction(s)
%            config.preprocess = 'rescale_diff';
%             config.preprocess_shift = [0 1];
%             x_states = featureNormailse(x_states,config);
%             y_states = featureNormailse(y_states,config);
%             z_states = featureNormailse(z_states,config);
%             
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
            
            % plot states as surface plot
            if config.plot_states
                subplot(2,2,2)
                plot(input_sequence)
                title('Input')
                
                subplot(2,2,3)
                plot(states)
                title('Scaled states')
                drawnow
                
                if ~config.parallel% surf plot of magnetisation direction
                    
                    subplot(2,2,4)
                    node_grid_size = sqrt(size(x_states,2));
                    [X,Y] = meshgrid(1:1:node_grid_size, 1:1:node_grid_size);
                    h = quiver3(X,Y,zeros(size(X)),reshape(x_states(1,1:end),node_grid_size,node_grid_size)...
                        ,reshape(y_states(1,1:end),node_grid_size,node_grid_size),...
                        reshape(z_states(1,1:end),node_grid_size,node_grid_size));
                    
                    set(gca, 'xDir', 'reverse',...
                        'camerapositionmode','manual','cameraposition',[1 1 1.25]);
                    
                    axis([1 node_grid_size 1 node_grid_size -1.25 1.25]);
                    zticks([-1 0 1])
                    zticklabels({min(min(states(:,1:end))) ,0, max(max(states(:,1:end)))})
                    grid off
                    
                    for t = 1:size(states,1)
                        %                         subplot(2,2,4)
                        %                         %newH = reshape(states{i}(t,:),node_grid_size,node_grid_size);
                        set(h,'UData',reshape(x_states(t,:),node_grid_size,node_grid_size));
                        set(h,'VData',reshape(y_states(t,:),node_grid_size,node_grid_size));
                        set(h,'WData',reshape(z_states(t,:),node_grid_size,node_grid_size));
                        
                        %                         subplot(2,2,1)
                        %                         Vq = interp2(reshape(z_states(t,1:end),node_grid_size,node_grid_size),5);
                        %                         imagesc(Vq)
                        %                         colormap(jet)
                        %                         drawnow
                    end
                end
            end
            
        else
            states = zeros(length(input_sequence),config.num_nodes(i)*length(config.read_mag_direction));
            fprintf('simulation failed\n')
            cmdout
        end
        
        if individual.core_indx(i) ~= 0 && config.num_reservoirs == 1
            % clean up files
            rmpath(file_path);
            
            try rmdir(file_path,'s');
                
            catch
               % warning('Error: Did not remove folder but contiued.\n')
            end
        end
    %end
    
end

%% change input signal to vampire simulator
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
    otherwise
        input = individual.input_scaling(1)*(individual.input_weights{1}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')';
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
    imagesc(rot90(reshape((max(abs(input))),node_grid_size,node_grid_size),3)')
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

%% change input file to vampire simulation
function changeInputFile(file_path, individual, input_sequence, indx, config)

% change parameters in input file
input_source=fopen(strcat(file_path,'default_input'),'r');
input_file=fopen(strcat(file_path,'input_file'),'w');

while ~feof(input_source)
    
    l=fgetl(input_source); % get line from base file, check if needs to be rewritten
    
    % crystal type
    if contains(l,'create:crystal-structure')
        switch(individual.material_type{1})
            case 'random_alloy'
                l = sprintf('create:crystal-structure = fcc'); % should always be fcc crystal
            otherwise
                l = sprintf('create:crystal-structure = %s', config.crystal_structure{1});
        end
    end
    % periodic boundaries
    if contains(l,'#create:periodic-boundaries-x') && individual.periodic_boundary(indx,1)
        l = sprintf('create:periodic-boundaries-x');
    end
    if contains(l,'#create:periodic-boundaries-y') && individual.periodic_boundary(indx,2)
        l = sprintf('create:periodic-boundaries-y');
    end
    if contains(l,'#create:periodic-boundaries-z') && individual.periodic_boundary(indx,3)
        l = sprintf('create:periodic-boundaries-z');
    end
    
    %core shell details
    if contains(l,'#create:sphere') && config.core_shell(indx)
        l = sprintf('create:sphere');
    end
    if contains(l,'#dimensions:particle-size') && (config.core_shell(indx) || ~strcmp(individual.material_shape{1},'film'))
        l = sprintf('dimensions:particle-size = %d !nm',individual.particle_size(indx));
    end
    
    %create shape
    if contains(l,'#create:sphere') && ~strcmp(individual.material_shape{1},'film')
        l = sprintf(strcat('create:',individual.material_shape{1}));
    end
    
    % system dimensions
    if contains(l,'dimensions:unit-cell-size')
        l = sprintf('dimensions:unit-cell-size = %.3f %s', config.unit_cell_size(1),config.unit_cell_units{1});
    end
    if contains(l,'dimensions:system-size-x')
        l = sprintf('dimensions:system-size-x = %d %s', individual.system_size(indx,1), config.size_units{1});
    end
    if contains(l,'dimensions:system-size-y')
        l = sprintf('dimensions:system-size-y = %d %s', individual.system_size(indx,2), config.size_units{2});
    end
    if contains(l,'dimensions:system-size-z')
        l = sprintf('dimensions:system-size-z = %.2f %s', individual.system_size(indx,3), config.size_units{3});
    end
    if contains(l,'cells:macro-cell-size')
        l = sprintf('cells:macro-cell-size = %d %s', individual.macro_cell_size(indx), config.macro_cell_units{1});
    end
    
    % sim details
    if contains(l,'sim:temperature')
        l = sprintf('sim:temperature = %d', individual.temperature(indx));
    end
    if contains(l,'sim:time-step')
        if contains(l,'sim:time-steps-increment')
            l = sprintf('sim:time-steps-increment = %d', individual.time_steps_increment);
        else
            l = sprintf('sim:time-step = %d %s', config.time_step,config.time_units);
        end
    end
    if contains(l,'sim:total-time-steps')
        config.total_time_steps = individual.time_steps_increment * size(input_sequence,1) + individual.time_steps_increment;
        l = sprintf('sim:total-time-steps = %d', config.total_time_steps);
    end
    % applied field
    if contains(l,'sim:applied-field-strength')
        l = sprintf('sim:applied-field-strength = %.2f !T' , individual.applied_field_strength(indx));
    end
    if contains(l,'sim:applied-field-unit-vector')
        l = sprintf('sim:applied-field-unit-vector = %s' , config.applied_field_unit_vector{1});
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
    
    fprintf(input_file,'%s \n',l);  % print line to file
end

fclose(input_source);
fclose(input_file);
end

%% change material file to vampire simulation
function changeMaterialFile(file_path, individual, indx, config)

% change parameters in material file
mat_file=fopen(strcat(file_path,'material_file'),'w');

switch(individual.material_type{indx})
    
    % TO DO!
    case 'multilayer'
        mat_source=fopen(strcat(file_path,'default_multi_material.txt'),'r');
        for m = 1:individual.num_materials(indx)
            while ~feof(mat_source)
                l=fgetl(mat_source); % get line from base file, check if needs to be rewritten
                if contains(l,strcat('material[',num2str(m),']:damping-constant='))
                    l = sprintf(strcat('material[',num2str(m),']:damping-constant=%d'), individual.damping(indx,m));
                end
                
                if contains(l,strcat('material[',num2str(m),']:exchange-matrix')) % need to apdapt this
                    for k = 1:individual.num_materials(indx)
                        if m == k
                            l = sprintf(strcat('material[',num2str(m),']:exchange-matrix[',num2str(k),']=%s'), individual.exchange(indx,m));
                        else
                            l = sprintf(strcat('material[',num2str(m),']:exchange-matrix[',num2str(k),']=%s'), individual.interfacial_exchange(indx));
                        end
                    end
                end
                
                if contains(l,strcat('material[',num2str(m),']:atomic-spin-moment'))
                    l = sprintf(strcat('material[',num2str(m),']:atomic-spin-moment=%d !muB'), individual.magmoment(indx,m));
                end
                if contains(l,strcat('material[',num2str(m),']:second-order-uniaxial-anisotropy-constant='))
                    l = sprintf(strcat('material[',num2str(m),']:second-order-uniaxial-anisotropy-constant=%s'), individual.anisotropy(indx,m));
                end
                if contains(l,strcat('material[',num2str(m),']:initial-spin-direction='))
                    l = sprintf(strcat('material[',num2str(m),']:initial-spin-direction=%s'), config.initial_spin_direction{indx});
                end
                if contains(l,strcat('material[',num2str(m),']:minimum-height'))
                    l = sprintf(strcat('material[',num2str(m),']:minimum-height=%.1f'), individual.minimum_height(indx,m));
                end
                if contains(l,strcat('material[',num2str(m),']:maximum-height'))
                    l = sprintf(strcat('material[',num2str(m),']:maximum-height=%.1f'), individual.maximum_height(indx,m));
                end
                fprintf(mat_file,'%s \n',l);  % print line to file
            end
        end
        
    case 'core_shell'
        
        mat_source=fopen(strcat(file_path,'default_core_shell.txt'),'r');
        for m = 1:individual.num_materials(indx)
            while ~feof(mat_source)
                l=fgetl(mat_source); % get line from base file, check if needs to be rewritten
                if contains(l,strcat('material[',num2str(m),']:damping-constant'))
                    l = sprintf(strcat('material[',num2str(m),']:damping-constant=%d'), individual.damping(indx,m));
                end
                if contains(l,strcat('material[',num2str(m),']:exchange-matrix')) % need to apdapt this
                    for k = 1:individual.num_materials(indx)
                        if m == k
                            l = sprintf(strcat('material[',num2str(m),']:exchange-matrix[',num2str(k),']=%s'), individual.exchange(indx,m));
                        else
                            l = sprintf(strcat('material[',num2str(m),']:exchange-matrix[',num2str(k),']=%s'), individual.interfacial_exchange(indx));
                        end
                    end
                end
                if contains(l,strcat('material[',num2str(m),']:atomic-spin-moment'))
                    l = sprintf(strcat('material[',num2str(m),']:atomic-spin-moment=%d !muB'), individual.magmoment(indx,m));
                end
                if contains(l,strcat('material[',num2str(m),']:second-order-uniaxial-anisotropy-constant='))
                    l = sprintf(strcat('material[',num2str(m),']:second-order-uniaxial-anisotropy-constant=%s'), individual.anisotropy(indx,m));
                end
                if contains(l,strcat('material[',num2str(m),']:initial-spin-direction'))
                    l = sprintf(strcat('material[',num2str(m),']:initial-spin-direction=%s'), config.initial_spin_direction{indx});
                end
                if contains(l,strcat('material[',num2str(m),']:core-shell-size'))
                    l = sprintf(strcat('material[',num2str(m),']:core-shell-size = %.1f'), individual.core_shell_size(indx,m));
                end
                fprintf(mat_file,'%s \n',l);  % print line to file
            end
        end
        
    case 'random_alloy'
        
        mat_source=fopen(strcat(file_path,'default_random_alloy.txt'),'r');
        for m = 1:individual.num_materials(indx)
            while ~feof(mat_source)
                l=fgetl(mat_source); % get line from base file, check if needs to be rewritten
                if contains(l,strcat('material[',num2str(m),']:damping-constant='))
                    l = sprintf(strcat('material[',num2str(m),']:damping-constant=%d'), individual.damping(indx,m));
                end
                
                if contains(l,strcat('material[1]:exchange-matrix[1]')) % need to apdapt this
                    l = sprintf(strcat('material[1]:exchange-matrix[1]=%s'), individual.exchange(indx,1));
                end
                if contains(l,strcat('material[1]:exchange-matrix[2]')) % need to apdapt this
                    l = sprintf(strcat('material[1]:exchange-matrix[2]=%s'), individual.interfacial_exchange(indx));
                end
                if contains(l,strcat('material[2]:exchange-matrix[1]')) % need to apdapt this
                    l = sprintf(strcat('material[2]:exchange-matrix[1]=%s'), individual.interfacial_exchange(indx));
                end
                if contains(l,strcat('material[2]:exchange-matrix[2]')) % need to apdapt this
                    l = sprintf(strcat('material[2]:exchange-matrix[2]=%s'), individual.exchange(indx,2));
                end
                
                if contains(l,strcat('material[',num2str(m),']:atomic-spin-moment'))
                    l = sprintf(strcat('material[',num2str(m),']:atomic-spin-moment=%d !muB'), individual.magmoment(indx,m));
                end
                if contains(l,strcat('material[',num2str(m),']:second-order-uniaxial-anisotropy-constant='))
                    l = sprintf(strcat('material[',num2str(m),']:second-order-uniaxial-anisotropy-constant=%s'), individual.anisotropy(indx,m));
                end
                if contains(l,strcat('material[',num2str(m),']:initial-spin-direction='))
                    l = sprintf(strcat('material[',num2str(m),']:initial-spin-direction=%s'), config.initial_spin_direction{indx});
                end
                if contains(l,'material[1]:alloy-fraction')
                    l = sprintf('material[1]:alloy-fraction[2]=%.2f', individual.alloy_fraction(indx));
                end
                fprintf(mat_file,'%s \n',l);  % print line to file
            end
        end
        
    otherwise
        
        mat_source=fopen(strcat(file_path,'default_material.txt'),'r');
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
            % finish...
            if contains(l,'material[1]:density')
                l = sprintf('material[1]:density=%.3f', individual.material_density(indx));
            end
            
            if contains(l,'material[1]:initial-spin-direction=')
                l = sprintf('material[1]:initial-spin-direction=%s', config.initial_spin_direction{1});
            end
            
            if contains(l,'material[1]:temperature-rescaling-exponent') && config.temperature_rescaling_exponent(2) > 0
                l = sprintf('material[1]:temperature-rescaling-exponent=%.3f', individual.temperature_rescaling_exponent(indx));
            end
            
            if contains(l,'material[1]:temperature-rescaling-curie-temperature') && config.temperature_rescaling_curie_temperature(2) > 0
                l = sprintf('material[1]:temperature-rescaling-curie-temperature=%d', individual.temperature_rescaling_curie_temperature(indx));
            end
            
            if contains(l,'material[1]:geometry-file') && ~isempty(config.geometry_file)
                l = sprintf('material[1]:geometry-file=%s',config.geometry_file);
                
                % find file and add to core
                try
                    geo_path = which('geometry_files.txt');
                catch
                    error('Error: Could not find the .geo file provided')
                end
                geo_path = geo_path(1:end-18);
                
                if config.evolve_geometry
                    % create new geo file
                    geo_file=fopen(strcat(file_path,config.geometry_file),'w');
                    
                    if config.evolve_poly
                        fprintf(geo_file,'%d \n',config.poly_num); % add number of lines
                        coord = individual.poly_coord;
                    else
                        % add new coordinates
                        fprintf(geo_file,'4 \n'); % add number of lines
                        coord = createRect(individual.geo_width,individual.geo_height); % get coords
%                         rectangle('Position',[0 0 individual.geo_width individual.geo_height])
%                         axis([0 1 0 1])
%                         drawnow
                    end
                    fprintf(geo_file,'%.4f %.4f \n',coord');
                    fclose(geo_file);
                else
                    % copy the known geo file to core directory
                    copyfile(strcat(geo_path,config.geometry_file),strcat(file_path,config.geometry_file));
                end
            end
            fprintf(mat_file,'%s \n',l);  % print line to file
        end
        
end

fclose(mat_source);
fclose(mat_file);

end