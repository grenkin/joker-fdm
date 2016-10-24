# Calculate the temperature and the radiative intensity in the 3-layer
#   complex heat transfer model
# Input data:
#   layer -- 1 x 3 vector of structures with fields: L, omega, n, Nc, A
#   boundary -- 1 x 2 vector of structures with fields: R, thetab, n
#   refraction -- logical value --
#     taking into account refraction and reflection at the interfaces
# Output data:
#   grid_info
#   phi, theta -- 1 x M cell matrices of 1 x K(j)+1 matrices
function [grid_info, theta, phi] = trilayer_calc (layer_data, boundary_data, refraction)
  for j = 1 : 3
    layer(j) = calc_layer(layer_data(j));
  endfor
  for j = 1 : 2
    boundary(j) = calc_boundary(boundary_data(j));
  endfor

  # u(1) -- theta
  # u(2) -- phi
  data.N = N = 2;
  data.M = M = 3;
  data.a(1, :) = [layer(:).C];
  data.a(2, :) = [layer(:).B];
  fm =  {@(x) x^4,  @(x) -x;
         @(x) -x^4, @(x) x};
  dfm = {@(x) 4*x^3,  @(x) -1;
         @(x) -4*x^3, @(x) 1};
  data.f = data.df = cell(N, M, N);
  for j = 1 : M
    for i = 1 : N
      for k = 1 : N
        data.f{i, j, k} = @(x) layer(j).K * fm{i, k}(x);
        data.df{i, j, k} = @(x) layer(j).K * dfm{i, k}(x);
      endfor
    endfor
  endfor
  data.b = [ Inf, Inf ; [boundary(:).tilde_gamma] ];
  data.v = [ [boundary(:).thetab] ; [boundary(:).thetab] .^ 4 ];
  if (refraction)
    for j = 1 : M - 1
      G_val(j) = calc_G(layer_data(j).n, layer_data(j + 1).n);
    endfor
  else
    G_val = Inf * ones(1, M - 1);
  endif
  data.G = [ Inf, Inf ; G_val ];
  L = [layer(:).L];
  K = 100 * ones(1, M);
  addpath("../../bvp1d")
  grid_info = get_grid_info(L, K);
  data.g = zeros(N, grid_info.nodes);
  guess = [ boundary(1).thetab + ...
    (boundary(2).thetab - boundary(1).thetab) * (0 : (grid_info.nodes - 1)) / ...
    (grid_info.nodes - 1) ; zeros(1, grid_info.nodes) ];
  tol = 1e-6;
  sol = solve_bvp1d(grid_info, data, guess, tol);
  sol_theta = sol(1, :);
  sol_phi = sol(2, :);

  theta = get_func_layers(grid_info, sol_theta);
  phi = get_func_layers(grid_info, sol_phi);
endfunction

function ret = calc_layer (l)
  alpha = 1 / (3 - l.A * l.omega);
  D = l.n ^ 2 * alpha;
  ret.B = D / l.L;
  ret.C = l.Nc * l.n ^ 2 / l.L;
  ret.K = l.L * l.n ^ 2 * (1 - l.omega);
  ret.L = l.L;
endfunction

function ret = calc_boundary (b)
  emiss = 1 - b.R;
  gam = gamma_func(emiss);
  ret.tilde_gamma = b.n ^ 2 * gam;
  ret.thetab = b.thetab;
endfunction

function ret = gamma_func (emiss)
  ret = emiss / ( 2 * (2 - emiss) );
endfunction
