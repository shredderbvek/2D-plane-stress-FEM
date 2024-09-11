%material properties
poisons_ratio = 0.3;
modulus_of_elasticity = 3 * 1e7;
%construction of D matrix
D_matrix =  zeros(3);
D_matrix(1,1) = (modulus_of_elasticity/(1 - poisons_ratio^2)) *1;
D_matrix(1,2) = (modulus_of_elasticity/(1 - poisons_ratio^2)) * poisons_ratio;
D_matrix(2,1) = (modulus_of_elasticity/(1 - poisons_ratio^2)) * poisons_ratio;
D_matrix(2,2) = (modulus_of_elasticity/(1 - poisons_ratio^2)) *1;
D_matrix(3,3) = (modulus_of_elasticity/(1 - poisons_ratio^2)) *(1 - poisons_ratio)/2;

%geometrical properties and mesh information

%global coordinates in the order 1, 2, 3, and 4
x = [0 ; 2 ; 2; 0];
y = [0, 0.5, 1, 1];
no_of_nodes = 4;
thickness = 1;

%dirilict boundary condition: displacement field

%first column denotes the node and the second column defines its displacement in meter
nodes_dirilict_values = [5 0; 6 0; 7 0; 8 0];
nodes_unknown = [3; 4; 5; 6];
displacement = zeros(no_of_nodes * 2, 1);

% Neumann boundary condition: traction

F = [0; 0; 0; 0; 0; -20; 0; -20];

%numerical integration using 4 point gauss quadrature
zeta_i = [-1/sqrt(3); 1/sqrt(3); -1/sqrt(3); 1/sqrt(3)];
eta_j = [-1/sqrt(3); -1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
weight_value = 1;

global_K = zeros(no_of_nodes *2);

for i = 1:4
    temp = weight_value^2 * k_function(zeta_i(i), eta_j(i), x, y, D_matrix);
    global_K = global_K + temp;
end

%computation of reduced system of equations

F_reduced = F(nodes_unknown);
K_reduced = global_K(nodes_unknown, nodes_unknown);

U_reduced = K_reduced\F_reduced;

global_K = global_K * thickness;

%assembly of total deformation
displacement(nodes_unknown) = U_reduced;

x_displacement = zeros(no_of_nodes,1);
y_displacement = zeros(no_of_nodes,1);

x_counter = 1;
y_counter = 1;

for i= 1:8
    if rem(i,2) == 0
        y_displacement(y_counter) = displacement(i);
        y_counter = y_counter + 1;
    else
        x_displacement(x_counter) = displacement(i);
        x_counter = x_counter + 1;
    end
end

%computes the magnitude of nodal displacements
u_total = sqrt(x_displacement.^2 + y_displacement.^2);

