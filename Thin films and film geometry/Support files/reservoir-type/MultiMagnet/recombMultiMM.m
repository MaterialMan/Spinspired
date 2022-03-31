%% Infection phase - need to update
function loser = recombMultiMM(winner,loser,config)

for layer = 1: config.num_layers
    
    W= winner.layer(layer).input_scaling(:);
    L = loser.layer(layer).input_scaling(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).input_scaling = reshape(L,size(loser.layer(layer).input_scaling));
    
    W= winner.layer(layer).leak_rate(:);
    L = loser.layer(layer).leak_rate(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).leak_rate = reshape(L,size(loser.layer(layer).leak_rate));
    
%     % params - W_scaling
%     W= winner.layer(layer).W_scaling(:);
%     L = loser.layer(layer).W_scaling(:);
%     pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
%     L(pos) = W(pos);
%     loser.layer(layer).W_scaling = reshape(L,size(loser.layer(layer).W_scaling));
%     
    %% magnet params
    W= winner.layer(layer).damping(:);
    L = loser.layer(layer).damping(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).damping = reshape(L,size(loser.layer(layer).damping));
    
    W= winner.layer(layer).anisotropy(:);
    L = loser.layer(layer).anisotropy(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).anisotropy = reshape(L,size(loser.layer(layer).anisotropy));
    
    W= winner.layer(layer).temperature(:);
    L = loser.layer(layer).temperature(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).temperature = reshape(L,size(loser.layer(layer).temperature));
    
    W= winner.layer(layer).exchange(:);
    L = loser.layer(layer).exchange(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).exchange = reshape(L,size(loser.layer(layer).exchange));
    
    W= winner.layer(layer).magmoment(:);
    L = loser.layer(layer).magmoment(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).magmoment = reshape(L,size(loser.layer(layer).magmoment));
    
    W= winner.layer(layer).applied_field_strength(:);
    L = loser.layer(layer).applied_field_strength(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).applied_field_strength = reshape(L,size(loser.layer(layer).applied_field_strength ));
    
    % layer thickness
    W= winner.layer(layer).thickness(:);
    L = loser.layer(layer).thickness(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).thickness = reshape(L,size(loser.layer(layer).thickness));
    
    % time_steps_increment
    W= winner.layer(layer).time_steps_increment(:);
    L = loser.layer(layer).time_steps_increment(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.layer(layer).time_steps_increment = reshape(L,size(loser.layer(layer).time_steps_increment));
    
    % material densities
    if config.evolve_material_density
        W= winner.layer(layer).material_density(:);
        L = loser.layer(layer).material_density(:);
        pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
        L(pos) = W(pos);
        loser.layer(layer).material_density = reshape(L,size(loser.layer(layer).material_density));
    end
    
    if config.evolve_geometry
        if config.evolve_poly
            W= winner.layer(layer).poly_coord(:);
            L = loser.layer(layer).poly_coord(:);
            pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
            L(pos) = W(pos);
            loser.layer(layer).poly_coord = reshape(L,size(loser.layer(layer).poly_coord));
        else
            W= winner.layer(layer).geo_width(:);
            L = loser.layer(layer).geo_width(:);
            pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
            L(pos) = W(pos);
            loser.layer(layer).geo_width = reshape(L,size(loser.layer(layer).geo_width));
            
            W= winner.layer(layer).geo_height(:);
            L = loser.layer(layer).geo_height(:);
            pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
            L(pos) = W(pos);
            loser.layer(layer).geo_height = reshape(L,size(loser.layer(layer).geo_height));
        end
    end
    
    
    %% inputs
    for i = 1:config.num_res_in_layer(layer)
        
        % input weights
        W= winner.layer(layer).input_weights{i}(:);
        L = loser.layer(layer).input_weights{i}(:);
        
%         % temporary for single input
%         ind_W= find(W); % find inputs in use
%         ind_L= find(L); % find inputs in use
%         
%         tmp_L = L(ind_L); %reset
%         L(ind_L) = 0;
%         posW = randperm(length(ind_W),ceil(config.rec_rate*length(ind_W)));
%         posL = randperm(length(ind_L),ceil(config.rec_rate*length(ind_L)));
%         tmp_L(posL) = W(ind_W(posW));
%         L(ind_W) = tmp_L;
        
        % take input
        L = W;
        
        % previous config
        % pos = randperm(length(L),ceil(config.rec_rate*length(L)));    %sum(rand(length(L),1) < config.rec_rate)
        % L(pos) = W(pos);
        
        loser.layer(layer).input_weights{i} = reshape(L,size(loser.layer(layer).input_weights{i}));
        
        % input widths
        if config.input_widths
        W= winner.layer(layer).input_widths{i}(:);
        L = loser.layer(layer).input_widths{i}(:);
        pos = randperm(length(L),ceil(config.rec_rate*length(L)));
        L(pos) = W(pos);
        loser.layer(layer).input_widths{i} = reshape(L,size(loser.layer(layer).input_widths{i}));       
        end
        
        % interfacial exchange
%         if config.random_alloy(i) || config.core_shell(i)
%             W= winner.layer(layer).interfacial_exchange(i);
%             L = loser.layer(layer).interfacial_exchange(i);
%             pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
%             L(pos) = W(pos);
%             loser.layer(layer).interfacial_exchange(i) = reshape(L,size(loser.layer(layer).interfacial_exchange(i)));
%         end
%         
%         if config.random_alloy(i)
%             % random alloy params
%             W= winner.layer(layer).alloy_fraction(i);
%             L = loser.layer(layer).alloy_fraction(i);
%             pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
%             L(pos) = W(pos);
%             loser.layer(layer).alloy_fraction(i) = reshape(L,size(loser.layer(layer).alloy_fraction(i)));
%         end
%         
%         if config.core_shell(i)
%             % core shell params
%             W= winner.layer(layer).shell_size(i,2);
%             L = loser.layer(layer).shell_size(i,2);
%             pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
%             L(pos) = W(pos);
%             loser.layer(layer).shell_size(i,2) = reshape(L,size(loser.layer(layer).shell_size(i,2)));
%         end
        
        % boundary params
        W= winner.layer(layer).periodic_boundary(i,logical(config.periodic_boundary));
        L = loser.layer(layer).periodic_boundary(i,logical(config.periodic_boundary));
        pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
        L(pos) = W(pos);
        loser.layer(layer).periodic_boundary(i,logical(config.periodic_boundary)) = reshape(L,size(loser.layer(layer).periodic_boundary(i,logical(config.periodic_boundary))));
        
        
        % layers
        loser.layer(layer).minimum_height(i,:) = [0 loser.layer(layer).thickness(i)];
        loser.layer(layer).maximum_height(i,:) = [loser.layer(layer).thickness(i) 1];
    end
    
    % for output weights
    if config.evolve_output_weights
        W= winner.output_weights(:);
        L = loser.output_weights(:);
        pos = randperm(length(L),ceil(config.rec_rate*length(L)));
        L(pos) = W(pos);
        loser.output_weights = reshape(L,size(loser.output_weights));
    end
    
end