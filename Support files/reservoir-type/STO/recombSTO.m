%% Infection phase - need to update
function loser = recombSTO(winner,loser,config)


for layer = 1:config.num_layers
    % reservoir parameters
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'input_scaling');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'leak_rate');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'W_scaling');
    
    %% magnet params
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'damping');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'anisotropy');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'temperature');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'exchange');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'magmoment');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'applied_field_strength');
    
    % time scales
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'interpolation_length');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'time_steps_increment');
    
    % spin torqe parameters
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'spin_transfer_relaxation_torque');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'spin_transfer_precession_torque');
    loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'spin_transfer_torque_asymmetry');
 
    % input weights

%     W= winner.layer(layer).W(:);
%     L = loser.layer(layer).W(:);
%     indices = find(loser.layer(layer).array_mask);
%     pos = randperm(length(indices),sum(rand(length(indices),1) < config.rec_rate));
%     L(indices(pos)) = W(indices(pos));
%     loser.layer(layer).W = reshape(L,size(loser.layer(layer).W));
    
    % remove weights
    for res_indx = 1:length(loser.layer(layer).nodes)
        loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'input_weights',res_indx);
        loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'W',res_indx);
        loser.layer(layer).W{res_indx}(~loser.layer(layer).array_mask{res_indx}) = 0;
    end
    
    % for output weights
    if config.evolve_output_weights
        loser.layer(layer) = recombFcn(winner.layer(layer),loser.layer(layer),config,'output_weights');
    end
end
end

function loser = recombFcn(winner,loser,config,parameter,indx)

if nargin == 5
    W = winner.(parameter){indx}(:);
    L = loser.(parameter){indx}(:);
else
    W = winner.(parameter)(:);
    L = loser.(parameter)(:);
end

pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);

if nargin == 5
    loser.(parameter){indx} = reshape(L,size(loser.(parameter){indx}));
else
    loser.(parameter) = reshape(L,size(loser.(parameter)));
end

end