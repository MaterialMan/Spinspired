function states = getLeakStates(states,individual,indx)
    leak_states = states(1,:);
    for n = 2:size(states,1)
        leak_states(n,:) = (1-individual.leak_rate(indx))*leak_states(n-1,:)+ individual.leak_rate(indx)*states(n,:);
    end
    states = leak_states;
end