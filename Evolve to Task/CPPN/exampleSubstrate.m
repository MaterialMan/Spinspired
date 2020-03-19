clear
close all

figure 

num_inputs = 2;
num_outputs = 2;
subrate_size = 6;
W = zeros(subrate_size.^2);

% define hidden nodes
[X_grid,Y_grid] = ndgrid(linspace(-1,1,subrate_size));

cnt = 1;
for i = 1:length(X_grid)
   for j = 1:length(Y_grid)
        hidden_node(cnt).X = X_grid(i,j);
        hidden_node(cnt).Y = Y_grid(i,j);
        hidden_node(cnt).Z = 0; 
        cnt = cnt+1;
   end
end
%scatter3([hidden_node.X],[hidden_node.Y],[hidden_node.Z])
hold on

% define input nodes
input_loc = linspace(-0.5,0.5,num_inputs);
cnt = 1;
for i = 1:length(input_loc)
        input_node(cnt).X = input_loc(i);
        input_node(cnt).Y = 0;
        input_node(cnt).Z = -1; 
        cnt = cnt+1;
end
%scatter3([input_node.X],[input_node.Y],[input_node.Z])

% define output nodes
output_loc = linspace(-0.5,0.5,num_outputs);
cnt = 1;
for i = 1:length(output_loc)
        output_node(cnt).X = output_loc(i);
        output_node(cnt).Y = 0;
        output_node(cnt).Z = 1; 
        cnt = cnt+1;
end
%scatter3([output_node.X],[output_node.Y],[output_node.Z])


% select node
node_to = 1;
node_from = 14;

% call fcn, e.g.
xy1 = [hidden_node(node_to).X hidden_node(node_to).Y hidden_node(node_to).Z];
xy2 = [hidden_node(node_from).X hidden_node(node_from).Y hidden_node(node_from).Z];
% call CPPN...

% return weight (example)
weight = rand;
% add to W matrix
W(node_to,node_from) = weight;
