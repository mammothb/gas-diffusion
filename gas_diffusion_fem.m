clear;
clf;

% Model parameters
para = Parameters();  % general parameters
gas_no = Gas('NO', para);  % create object for NO gas
gas_o2 = Gas('O2', para);  % create object for O2 gas

cfl = 3;

elem_per_compt = 5;  % Number of elements per compartment
num_compt = 5;  % Number of compartments
nodes_per_elem = 2;  % Number of (local) nodes per element

num_elems = num_compt * elem_per_compt;  % Total number of elements
num_nodes = num_elems + 1;  % Total number of nodes

node_pos = zeros(num_nodes, 1);  % position of node n
elem_node = zeros(num_nodes, 2);  % global node coordinate of element e
bc = zeros(num_nodes, 2);  % boundary conditions at node n

% Different space steps in each compartment
h = zeros(num_compt, 1);
h(1) = (para.int_r - cfl) / elem_per_compt;
h(2) = cfl / elem_per_compt;
h(3) = para.len_EC / elem_per_compt;
h(4) = para.len_VW / elem_per_compt;
h(5) = para.len_T / elem_per_compt;
disp(para.R);

% Create mesh
compt_ind = 1;
for ii = 1 : num_nodes
  if ii > 1
    node_pos(ii) = node_pos(ii - 1) + h(compt_ind);
    if mod(ii - 1, elem_per_compt) == 0
      compt_ind = compt_ind + 1;
    end
  end
  bc(ii, :) = [2, 0];
end

plot(node_pos, zeros(num_nodes, 1), '.');

% Generate elements
for ii = 1 : num_elems
  elem_node(ii, :) = [ii, ii + 1];
end
