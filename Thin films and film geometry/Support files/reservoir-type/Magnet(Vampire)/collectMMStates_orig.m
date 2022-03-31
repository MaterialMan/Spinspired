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
        
        % change/update source/input files
        if i > 1
            changeSourceFile(file_path, individual, input_sequence, i, states{i-1},config)
        else
            changeSourceFile(file_path, individual, input_sequence, i, states{i},config)
        end
        
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
        %command = strcat('cd "', file_path,'" && ./vampire-parallel');
        
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
            x_states = tmp_states(1:size(input_sequence,1),:);
            y_states = tmp_states(size(input_sequence,1)+1:size(input_sequence,1)*2,:);
            z_states = tmp_states(size(input_sequence,1)*2 +1:end,:);
            
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
                
                subplot(2,2,2)
                plot(input_sequence(config.wash_out+1:end,:))
                title('Input')
                
                subplot(2,2,3)
                plot(states{i}(config.wash_out+1:end,:))
                title('Scaled states')
                drawnow
                
                if config.parallel% surf plot of magnetisation direction
                    
                    subplot(2,2,4)
                    node_grid_size = sqrt(size(x_states,2));
                    [X,Y] = meshgrid(1:1:node_grid_size, 1:1:node_grid_size);
                    Z=zeros(size(X));
                    
                    f = getframe(gca);
                    input_indx = find(individual.input_weights{1}(:,1));
                    bias_indx = find(individual.input_weights{1}(:,2));
                    other_indx = 1:size(states{i},2)/3;
                    other_indx(input_indx) = 0;
                    other_indx(bias_indx) = 0;
                    other_indx(other_indx==0) = [];

%                     h = surf(X,Y,reshape(z_states(1,:),node_grid_size,node_grid_size));
%                     colormap(gca,'bone'); %
%                     %set(gca,'visible','off')
%                     set(gca,'XColor', 'none','YColor','none')
%                     shading interp
                    for t = 2:size(states{i},1)
                        
                      %set(h,'zdata',reshape(z_states(1,:),node_grid_size,node_grid_size),'facealpha',0.65);
                        
                        quiver3(X,Y,zeros(size(X)),zeros(size(X))...
                            ,zeros(size(X)),...
                            zeros(size(X)));
                        
                        hold on
                        quiver3(X(other_indx),Y(other_indx),Z(other_indx),x_states(t,other_indx)...
                            ,y_states(t,other_indx),...
                            z_states(t,other_indx),0.5,'Color','b');
                        
                        
                        quiver3(X(input_indx),Y(input_indx),Z(input_indx),x_states(t,input_indx)'...
                            ,y_states(t,input_indx)',...
                            z_states(t,input_indx)',0.5,'Color','r');

                        quiver3(X(bias_indx),Y(bias_indx),Z(bias_indx),x_states(t,bias_indx)'...
                            ,y_states(t,bias_indx)',...
                            z_states(t,bias_indx)',0.5,'Color','g');

                        set(gca, 'xDir', 'reverse',...
                            'camerapositionmode','manual','cameraposition',[1 1 1.25]);
                        
                        axis([1 node_grid_size 1 node_grid_size -1.25 1.25]);
                        zticks([-1 0 1])
                        zticklabels({min(min(states{i}(:,1:end))) ,0, max(max(states{i}(:,1:end)))})
                        grid off
                        hold off
                        drawnow
                        f(t) = getframe(gca);
                    end
                end
            end
        else
            states{i} = zeros(length(input_sequence),config.num_nodes(i));
            fprintf('simulation failed\n')
            cmdout
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
        
%     case 'pipeline_IA'
%         if indx < 2
%             input = individual.input_scaling(indx)*(individual.input_weights{indx}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')';
%         else
%             input = previous_states(:,1:end-individual.n_input_units*(size(previous_states,2) > individual.nodes(indx))) * (individual.input_scaling(indx)*individual.input_weights{indx})...
%                 + individual.W{indx-1,indx}*individual.W_scaling(indx-1,indx)*previous_states;
%             
%            % input = individual.input_scaling(indx)*(individual.input_weights{indx}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')'+...
%            % previous_states(:,1:end-individual.n_input_units*(size(previous_states,2) > individual.nodes(indx))) * (individual.W_scaling(indx-1,indx)*individual.W{indx-1,indx});
%         
%         end
    case 'pipeline'
        if indx < 2
            input = individual.input_scaling(indx)*(individual.input_weights{indx}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')';
        else
           % input =  [previous_states(:,1:end-individual.n_input_units*(size(previous_states,2) > individual.nodes(indx))) input_sequence repmat(individual.bias_node,size(input_sequence,1),1)] * (individual.W_scaling(indx-1,indx)*individual.W{indx-1,indx});
           if config.add_pipeline_input
               input =  [[zeros(1, size(previous_states,2)); previous_states(1:end-1,:)] input_sequence repmat(individual.bias_node,size(input_sequence,1),1)] * (individual.W_scaling(indx-1,indx)*individual.W{indx-1,indx}');
           else
               input =  [[zeros(1, size(previous_states,2)); previous_states(1:end-1,:)] repmat(individual.bias_node,size(input_sequence,1),1)] * (individual.W_scaling(indx-1,indx)*individual.W{indx-1,indx}');
           end
        end
        
    otherwise
         input = individual.input_scaling(indx)*(individual.input_weights{indx}*[input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]')';

end

% change input widths
if config.input_widths > 0
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
        switch(individual.material_type{indx})
            case 'random_alloy'
                l = sprintf('create:crystal-structure = fcc'); % should always be fcc crystal
            otherwise
                l = sprintf('create:crystal-structure = %s', config.crystal_structure{indx});
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
    if contains(l,'#dimensions:particle-size') && (config.core_shell(indx) || ~strcmp(individual.material_shape{indx},'film'))
        l = sprintf('dimensions:particle-size = %d !nm',individual.particle_size(indx));
    end
    
    %create shape
    if contains(l,'#create:sphere') && ~strcmp(individual.material_shape{indx},'film')
        l = sprintf(strcat('create:',individual.material_shape{indx}));
    end
    
    % system dimensions
    if contains(l,'dimensions:unit-cell-size')
        l = sprintf('dimensions:unit-cell-size = %.3f %s', config.unit_cell_size(indx),config.unit_cell_units{indx});
    end
    if contains(l,'dimensions:system-size-x')
        l = sprintf('dimensions:system-size-x = %d %s', individual.system_size(indx,1), config.size_units{indx,1});
    end
    if contains(l,'dimensions:system-size-y')
        l = sprintf('dimensions:system-size-y = %d %s', individual.system_size(indx,2), config.size_units{indx,2});
    end
    if contains(l,'dimensions:system-size-z')
        l = sprintf('dimensions:system-size-z = %d %s', individual.system_size(indx,3), config.size_units{indx,3}); %was.2f
    end
    if contains(l,'cells:macro-cell-size')
        l = sprintf('cells:macro-cell-size = %d %s', config.macro_cell_size(indx), config.macro_cell_units{indx});
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
        config.total_time_steps = individual.time_steps_increment(indx) * size(input_sequence,1) + individual.time_steps_increment(indx);
        l = sprintf('sim:total-time-steps = %d', config.total_time_steps);
    end
    % applied field
    if contains(l,'sim:applied-field-strength') 
        l = sprintf('sim:applied-field-strength = %.2f !T' , individual.applied_field_strength(indx));
    end
     if contains(l,'sim:applied-field-unit-vector') 
        l = sprintf('sim:applied-field-unit-vector = %s' , config.applied_field_unit_vector{indx});
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
            
            if contains(l,'material[1]:damping-constant')
                l = sprintf('material[1]:damping-constant=%.3f', individual.damping(indx));
            end
            if contains(l,'material[1]:exchange-matrix[1]')
                l = sprintf('material[1]:exchange-matrix[1]=%s', individual.exchange(indx));
            end
            if contains(l,'material[1]:atomic-spin-moment')
                l = sprintf('material[1]:atomic-spin-moment=%.3f !muB', individual.magmoment(indx));
            end
            if contains(l,'material[1]:second-order-uniaxial-anisotropy-constant')
                l = sprintf('material[1]:second-order-uniaxial-anisotropy-constant=%s', individual.anisotropy(indx));
            end
            % finish...
            if contains(l,'material[1]:density')
                l = sprintf('material[1]:density=%.3f', individual.material_density(indx));
            end
            
            if contains(l,'material[1]:initial-spin-direction')
                l = sprintf('material[1]:initial-spin-direction=%s', config.initial_spin_direction{indx});
            end
            
%             if contains(l,'material[1]:temperature-rescaling-exponent') && config.temperature_rescaling_exponent(2) > 0
%                 l = sprintf('material[1]:temperature-rescaling-exponent=%.3f', individual.temperature_rescaling_exponent(indx));
%             end
%             
%             if contains(l,'material[1]:temperature-rescaling-curie-temperature') && config.temperature_rescaling_curie_temperature(2) > 0
%                 l = sprintf('material[1]:temperature-rescaling-curie-temperature=%d', individual.temperature_rescaling_curie_temperature(indx));
%             end
            
            fprintf(mat_file,'%s \n',l);  % print line to file
        end
        
end

fclose(mat_source);
fclose(mat_file);

end