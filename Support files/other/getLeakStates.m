function states = getLeakStates(states,individual)
    leak_states = states(1,:);
    for n = 2:size(states,1)
        leak_states(n,:) = (1-individual.leak_rate)*leak_states(n-1,:)+ individual.leak_rate*states(n,:);
    end
    states = leak_states;
end