function states = getLeakStates(states,individual,config)

for n = 2:size(states,1)    
    if config.leak_delay
        if n > max(individual.max_update_cycle)
            ind = sub2ind(size(states),n-individual.update_cycle,1:size(states,2));
            tmp_states = states(ind);
            states(n,:) = (1-individual.leak_rate).*tmp_states + individual.leak_rate.*states(n,:);
        end
    else
        states(n,:) = (1-individual.leak_rate).*states(n-1,:) + individual.leak_rate.*states(n,:);
    end
end


