%{
Computes B_matrix for bilinear quadratic element at the given reference
coordinate
%}
function [result] = compute_B_matrix(zeta_val, eta_val, x, y)
local_coords = [-1, -1;1, -1;1, 1;-1, 1 ];
%construction of shape functions for bilinear quadrilatera element using
%isoparametric formulation
syms zeta eta
shape_functions = sym(zeros(1,4));  % Initialize symbolic array
for i = 1:4
    shape_functions(i) = 0.25*(1 + zeta * local_coords(i,1))* (1 + eta * local_coords(i,2));
end
%partial derivatives of shape functions
diff_N = jacobian(shape_functions, [zeta, eta]);
% mapping of physical domain with natural domain
syms x1 x2 x3 x4 y1 y2 y3 y4
real_2_ref_map = sym([x1 * shape_functions(1) + x2 * shape_functions(2)+ x3 * shape_functions(3)+ x4 * shape_functions(4); y1 * shape_functions(1) + y2 * shape_functions(2)+ y3 * shape_functions(3)+ y4 * shape_functions(4)]);
%computation of jacobian matrix
J_matrix = jacobian(real_2_ref_map, [zeta, eta]);
J_matrix = J_matrix';
%computation of inverse of jacobian matrix
J_inverse = inv(J_matrix);
%computation of A matrix
A_matrix = sym(zeros(4));
A_matrix(1, 1) = J_inverse(1,1);
A_matrix(1, 2) = J_inverse(1,2);
A_matrix(2, 1) = J_inverse(2,1);
A_matrix(2, 2) = J_inverse(2,2);
A_matrix(3, 3) = J_inverse(1,1);
A_matrix(3, 4) = J_inverse(1,2);
A_matrix(4, 3) = J_inverse(2,1);
A_matrix(4, 4) = J_inverse(2,2);
%computation of G matrix
G_matrix = sym(zeros(4,8));

for i = 1:4
    if i == 1
        k = 1;
        for j = 1:8
            if rem(j,2) == 1
                G_matrix(i, j) = diff_N(k,1);
                k = k + 1;
            end
        end
    end
    
    if i == 2
        k = 1;
        for j = 1:8
            if rem(j,2) == 1
                G_matrix(i, j) = diff_N(k,2);
                k = k + 1;
            end
        end
    end
    
    if i == 3
        k = 1;
        for j = 1:8
            if rem(j,2) == 0
                G_matrix(i, j) = diff_N(k,1);
                k = k + 1;
            end
        end
    end
    
    if i == 4
        k = 1;
        for j = 1:8
            if rem(j,2) == 0
                G_matrix(i, j) = diff_N(k,2);
                k = k + 1;
            end
        end
    end
end
% initialization of C matrix
C_matrix = [1, 0, 0, 0; 0, 0, 0, 1; 0, 1, 1, 0];
% computation of B matrix
B_matrix =  C_matrix * A_matrix * G_matrix;
result = double(subs(B_matrix, {zeta, eta, x1, x2, x3, x4, y1, y2, y3, y4}, {zeta_val, eta_val, x(1), x(2), x(3), x(4), y(1), y(2), y(3), y(4)}));
end