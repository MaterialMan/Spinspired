%% called from getDataSetInfo
% Create graph structure G to use

function [config,new_num_nodes,G] = getShape(config,nodes,type)

for graph_indx = 1:length(nodes)
    
    if contains(config.graph_type,'Lattice')
        if (sqrt(nodes(graph_indx)) ~= round(sqrt(nodes(graph_indx))))
            error('\n Number of nodes needs to be a square number. \n')
        else
            nodes(graph_indx) = sqrt(nodes(graph_indx));
        end
    end
    
    num_nodes= nodes(graph_indx);
    
    if length(config.graph_type) > 1
        graph_type = config.graph_type{graph_indx};
    else
        graph_type = config.graph_type{1};
    end
    
    switch(graph_type)
        
        case 'Bucky'
            G{graph_indx} = graph(bucky);
            config.plot_3d = 1;      % plot graph in 3D.
            
        case 'L'
            A = delsq(numgrid('L',num_nodes +2));
            G{graph_indx} = graph(A,'omitselfloops');
            config.plot_3d = 0;    % plot graph in 3D.
            
        case 'Hypercube'
            A = hypercube(num_nodes);
            G{graph_indx} = graph(A);
            config.plot_3d = 1;    % plot graph in 3D.
            
        case 'Torus'
            config.rule_type = 'Moores';
            config.torus_rings = config.num_nodes;
            G{graph_indx} = torusGraph(num_nodes,config.self_loop(graph_indx),config);
            config.plot_3d = 1;    % plot graph in 3D.
            
        case 'Barbell'
            load barbellgraph.mat
            G{graph_indx} = graph(A,'omitselfloops');
            config.plot_3d = 0;    % plot graph in 3D.
            
        case {'basicLattice','partialLattice','fullLattice','basicCube','partialCube','fullCube','ensembleLattice', 'ensembleCube','ensembleShape'}
            G{graph_indx} = createLattice(num_nodes,graph_type,config.self_loop,1);
            config.plot_3d = 0;    % plot graph in 3D.
            
        case 'Ring'
            config.rule_type = 0;
            config.torus_rings = 1;
            G{graph_indx} = torusGraph(num_nodes,config.self_loop(graph_indx),config);
            config.plot_3d = 0;    % plot graph in 3D.
            
        otherwise
            error('Requires a substrate shape. Check graph type.')
    end
    
    %config.G{graph_indx} = G;
    new_num_nodes(graph_indx) = size(G{graph_indx}.Nodes,1);
    
    if contains(type,'diagraph')
        G{graph_indx} = digraph(adjacency(G{graph_indx}));  %convert to directed graph  
    end
end
