
rng(1,'twister')

indvidual.lattice_size = 200;
config.time_steps = 1000;
config.type = 'CML';
config.global = 1;

input_sequence = rand(time_steps,1);

%for i = 1:1000 
    
    indvidual.Win = zeros(indvidual.lattice_size,1);
    
    if ~config.global
        indvidual.W = rand(indvidual.lattice_size,1);
    else
        indvidual.W = ones(indvidual.lattice_size,1)*0.3;
    end
    
    indvidual.r = 1.75;
    indvidual.x0 = zeros(1,indvidual.lattice_size);

    CML(indvidual, input_sequence,config);
%end