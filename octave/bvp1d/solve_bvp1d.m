# Solve a nonlinear boundary-value problem in 1D

# grid_info -- information on the grid returned by get_grid_info
# data -- structure with fields:
#   N -- number of equations, M -- number of subdomains
#   a -- diffusion coefficients
#   f -- nonlinear terms functions, df -- derivatives of f,
#   g -- right-hand sides
#   b -- coefficients in BCs, v -- functions values in BCs
#   G -- coefficients in conjugation conditions
# guess -- initial guess
# tol -- tolerance in Newton's method

# Returned value (sol) -- values of N grid functions at grid points
# sol -- N x grid_info.nodes matrix

# a -- N x M matrix, a(i, j) > 0
# f, df -- N x M x N cell matrices of function handles
# g -- N x grid_info.nodes matrix
# b, v -- N x 2 matrices, b(i, j) >= 0 or b(i, j) == Inf
# G -- N x (M - 1) matrix, G(i, j) >= 0 or G(i, j) == Inf
# guess -- N x grid_info.nodes matrix

function sol = solve_bvp1d (grid_info, data, guess, tol)
  [N, M] = deal(data.N, data.M);
  unknowns = N * grid_info.nodes;
  A = sparse(unknowns, unknowns);
  rhs = zeros(unknowns, 1);
  # Apply Newton's method for solving the nonlinear system
  iter = 0;
  sol = guess;
  do
    printf("%d ", ++iter);
    sol_old = sol;
    # Build the matrix and the right-hand side of the linear system
    for i = 1 : N
      for j = 1 : M
        for n = 0 : grid_info.K(j)
          ind = nindex(grid_info, i, j, n);
          [A(ind, :), rhs(ind)] = ...
            get_fd_equation(grid_info, data, sol_old, i, j, n);
        endfor
      endfor
    endfor
    # Solve the linear system
    x = A \ rhs;
    sol = transpose(reshape(x, grid_info.nodes, N));
  until (norm(sol - sol_old) < tol)
  printf("\n");
endfunction

# Get FD equation in a given grid node
function [coeff, rhs] = get_fd_equation (grid_info, data, guess, i, j, n)
  M = data.M;
  if (j == 1 && n == 0 || j == M && n == grid_info.K(j))
    # boundary node
    [coeff, rhs] = get_eq_boundary(grid_info, data, guess, i, j, n);
  elseif (n == grid_info.K(j))
    # right-most node of the j-th domain
    [coeff, rhs] = get_eq_interface_right(grid_info, data, guess, i, j, n);
  elseif (n == 0)
    # left-most node of the j-th domain
    [coeff, rhs] = get_eq_interface_left(grid_info, data, guess, i, j, n);
  else
    # internal node
    [coeff, rhs] = get_eq_internal(grid_info, data, guess, i, j, n);
  endif
endfunction

# Get FD equation in an internal node
# 0 < n < data.K(j)
function [coeff, rhs] = get_eq_internal (grid_info, data, guess, i, j, n)
  N = data.N;
  # [coeff1, rhs1] represents the differential operator
  # [coeff2, rhs2] represents the nonlinear operator
  unknowns = N * grid_info.nodes;
  coeff1 = coeff2 = sparse(1, unknowns);
  ind = arrayfun(@(nn) nindex(grid_info, i, j, nn), [n - 1, n, n + 1]);  # stencil
  coeff1(ind) = - data.a(i, j) / grid_info.h(j) ^ 2 * [1, -2, 1];
  rhs1 = 0;
  [coeff2, rhs2] = get_nonlinear_term(grid_info, data, guess, i, j, n);
  [coeff, rhs] = deal(coeff1 + coeff2, rhs1 + rhs2);
endfunction

# Get FD equation in a boundary node
# (j == 1 and n == 0) or (j == M and n == grid_info.K(M))
function [coeff, rhs] = get_eq_boundary (grid_info, data, guess, i, j, n)
  [N, M] = deal(data.N, data.M);
  if (n == 0)
    # left boundary node
    n1 = n + 1;
    bval = data.b(i, 1);
    vval = data.v(i, 1);
  else
    # right boundary node
    # n == grid_info.K(M)
    n1 = n - 1;
    bval = data.b(i, 2);
    vval = data.v(i, 2);
  endif

  # [coeff1, rhs1] represents the differential operator
  # [coeff2, rhs2] represents the nonlinear operator
  unknowns = N * grid_info.nodes;
  coeff1 = coeff2 = sparse(1, unknowns);
  if (bval == Inf)
    # Dirichlet BC
    ind = nindex(grid_info, i, j, n);
    coeff1(ind) = 1;
    rhs1 = vval;
  else
    # Neumann or Robin BC
    ind = arrayfun(@(nn) nindex(grid_info, i, j, nn), [n, n1]);  # stencil
    coeff1(ind) = data.a(i, j) / grid_info.h(j) * [1, -1] + bval * [1, 0];
    rhs1 = bval * vval;
  endif
  [coeff2, rhs2] = get_nonlinear_term(grid_info, data, guess, i, j, n);
  coeff2 *= grid_info.h(j) / 2;
  rhs2 *= grid_info.h(j) / 2;
  [coeff, rhs] = deal(coeff1 + coeff2, rhs1 + rhs2);
endfunction

# Get FD equation in a right-most node of the j-th domain
# (interface node from the left between domains j and j+1)
# 1 <= j < M, n == grid_info.K(j)
function [coeff, rhs] = get_eq_interface_right (grid_info, data, guess, i, j, n)
  N = data.N;
  # [coeff1, rhs1] represents the differential operator
  # [coeff2, rhs2] represents the nonlinear operator
  unknowns = N * grid_info.nodes;
  coeff1 = coeff2 = sparse(1, unknowns);
  get_ind = @(jj, nn) nindex(grid_info, i, jj, nn);
  Gval = data.G(i, j);
  if (Gval == Inf)
    ind = [get_ind(j, n), get_ind(j, n - 1), get_ind(j + 1, 0), get_ind(j + 1, 1)];  # stencil
    coeff1(ind) = data.a(i, j) / grid_info.h(j) * [1, -1, 0, 0] + ...
      data.a(i, j + 1) / grid_info.h(j + 1) * [0, 0, 1, -1];
    rhs1 = 0;
    [coeff2a, rhs2a] = get_nonlinear_term(grid_info, data, guess, i, j, n);
    [coeff2b, rhs2b] = get_nonlinear_term(grid_info, data, guess, i, j + 1, 0);
    coeff2 = coeff2a * grid_info.h(j) / 2 + coeff2b * grid_info.h(j + 1) / 2;
    rhs2 = rhs2a * grid_info.h(j) / 2 + rhs2b * grid_info.h(j + 1) / 2;
  else
    ind = [get_ind(j, n), get_ind(j, n - 1), get_ind(j + 1, 0)];  # stencil
    coeff1(ind) = data.a(i, j) / grid_info.h(j) * [1, -1, 0] + Gval * [1, 0, -1];
    rhs1 = 0;
    [coeff2, rhs2] = get_nonlinear_term(grid_info, data, guess, i, j, n);
    coeff2 *= grid_info.h(j) / 2;
    rhs2 *= grid_info.h(j) / 2;
  endif
  [coeff, rhs] = deal(coeff1 + coeff2, rhs1 + rhs2);
endfunction

# Get FD equation in a left-most node of the j-th domain
# (interface node from the right between domains j-1 and j)
# 1 < j <= M, n == 0
function [coeff, rhs] = get_eq_interface_left (grid_info, data, guess, i, j, n)
  N = data.N;
  # [coeff1, rhs1] represents the differential operator
  # [coeff2, rhs2] represents the nonlinear operator
  unknowns = N * grid_info.nodes;
  coeff1 = coeff2 = sparse(1, unknowns);
  get_ind = @(jj, nn) nindex(grid_info, i, jj, nn);
  Gval = data.G(i, j - 1);
  if (Gval == Inf)
    ind = [get_ind(j, n), get_ind(j - 1, grid_info.K(j - 1))];  # stencil
    coeff1(ind) = [1, -1];
    rhs1 = 0;
    # coeff2 == 0
    rhs2 = 0;
  else
    ind = [get_ind(j, n), get_ind(j, n + 1), get_ind(j - 1, grid_info.K(j - 1))];  # stencil
    coeff1(ind) = data.a(i, j) / grid_info.h(j) * [1, -1, 0] + Gval * [1, 0, -1];
    rhs1 = 0;
    [coeff2, rhs2] = get_nonlinear_term(grid_info, data, guess, i, j, n);
    coeff2 *= grid_info.h(j) / 2;
    rhs2 *= grid_info.h(j) / 2;
  endif
  [coeff, rhs] = deal(coeff1 + coeff2, rhs1 + rhs2);
endfunction

# Get the linearized nonlinear term \sum_k( f_{ijk}(u_{kjn}) ) - g_{ijn}
function [coeff, rhs] = get_nonlinear_term (grid_info, data, guess, i, j, n)
  N = data.N;
  unknowns = N * grid_info.nodes;
  coeff = sparse(1, unknowns);
  gind = gindex(grid_info, j, n);
  for k = 1 : N
    ind = nindex(grid_info, k, j, n);
    coeff(ind) = data.df{i, j, k}( guess(k, gind) );
  endfor
  rhs = - data.g(i, gind);
  for k = 1 : N
    rhs += data.df{i, j, k}( guess(k, gind) ) * guess(k, gind) ...
    - data.f{i, j, k}( guess(k, gind) );
  endfor
endfunction

# Index of a grid node in a vector of unknowns
# i -- function index
# j -- subdomain index
# n -- grid node index within the subdomain
# 1 <= i <= N, 1 <= j <= M, 0 <= n <= K(j)
# N -- number of equations, M -- number of subdomains
# K(j) -- number of grid nodes in the j-th subdomain
function ret = nindex (grid_info, i, j, n)
  ret = (i - 1) * grid_info.nodes + gindex(grid_info, j, n);
endfunction
