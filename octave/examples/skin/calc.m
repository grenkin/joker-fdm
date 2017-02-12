# Radiation frequences
global FREQ_UV = 310  FREQ_VIS = 514  FREQ_IR = 800  # nanometers
# Taking into account the Planck function term (theta^4)
global PLANCK_OFF = 0  PLANCK_ON = 1
# Taking into account refraction and reflection at the interfaces
global REFR_OFF = 0  REFR_ON = 1
# Taking into account administration of nanoparticles
global NANO_OFF = 0  NANO_ON = 1

addpath("../joker-fdm/octave/bvp1d")

# Calculate the temperature and the radiative intensity in the skin model
# Returned values: xgrid, temp, intens, Q -- 1 x M cell matrices of 1 x grid_info.nodes matrices
function [grid_info, temp, intens, Q] = skin_calc (freq, planck, refr, nano)
  global FREQ_UV FREQ_VIS FREQ_IR PLANCK_OFF PLANCK_ON REFR_OFF REFR_ON NANO_OFF NANO_ON
  switch (freq)
    case FREQ_UV
      n = [1.5503, 1.53, 1.4, 1.4];
      mu_a = 100 * [3000, 600, 300, 8.7];
      mu_s = 100 * [4560, 2400, 1400, 583];
      g = [0.7, 0.9, 0.71, 0.71];    
    case FREQ_VIS
      n = [1.5426, 1.53, 1.4, 1.4];
      mu_a = 100 * [60, 60, 44, 2.2];
      mu_s = 100 * [2200, 1560, 600, 250];
      g = [0.7, 0.9, 0.77, 0.77];
    case FREQ_IR
      n = [1.5408, 1.53, 1.4, 1.4];
      mu_a = 100 * [3, 3, 40, 1.7];
      mu_s = 100 * [500, 420, 420, 175];
      g = [0.76, 0.9, 0.85, 0.85];
  endswitch
  k_cond = [0.266, 0.266, 0.498, 0.498];
  switch (nano)
    case NANO_ON
      M = 4;
      L = 1e-6 * [3, 17, 100, 500];
      K = [50, 100, 150, 200];
    case NANO_OFF
      M = 3;
      L = 1e-6 * [20, 100, 500];
      K = [150, 150, 200];
      n = prepad(n, M);
      mu_a = prepad(mu_a, M);
      mu_s = prepad(mu_s, M);
      g = prepad(g, M);
      k_cond = prepad(k_cond, M);
  endswitch

  StBol = 5.67e-8;
  Kelvin0 = 273.15;
  # Tmax = 1
  # phi(x) = ( pi / (StBol * n^2) ) * 1 / 2 * \int_{-1}^1 I(x,\nu) d\nu
  I_norm_coeff = pi / StBol;  # normalization coefficient

  alpha = 1 ./ (3 .* (mu_a + (1 .- g) .* mu_s));
  sigma = 4 * StBol * mu_a .* n .^ 2;

  # parameters of perfusion
  omega = 0.24;
  c_b = 4186;
  rho_b = 1000;
  theta_b = 37 + Kelvin0;

  # parameters of BC at the left point
  A = 100;
  theta_ext = 20 + Kelvin0;
  phi0 = 1000 * I_norm_coeff;
  Phi0 = 0 * I_norm_coeff;

  switch (refr)
    case REFR_ON
      G = zeros(1, M - 1);
      for i = 1 : M - 1
        G(i) = calc_G(n(i), n(i + 1));
      endfor      
    case REFR_OFF
      G = Inf * ones(1, M - 1);
  endswitch

  BC_coeff1 = calc_BC_coeff1(n(1));
  BC_coeff2 = calc_BC_coeff2(n(1));

  # u(1) -- phi
  # u(2) -- theta
  data.N = N = 2;
  data.M = M;
  data.a(1, :) = n .* alpha;
  data.a(2, :) = k_cond;
  data.b = [BC_coeff1, n(M) ^ 2 / 2 ; A, Inf];
  switch (planck)
    case PLANCK_ON
      blood_emiss = theta_b ^ 4;
    case PLANCK_OFF
      blood_emiss = 0;
  endswitch
  data.v = [phi0 + BC_coeff2 * Phi0 / BC_coeff1, blood_emiss ; theta_ext, theta_b];
  data.G = [ G ; Inf * ones(1, M - 1) ];

  grid_info = get_grid_info(L, K);

  switch (planck)
    case PLANCK_ON
      fm =  {@(x) x,  @(x) -x^4;
             @(x) -x, @(x) x^4};
      dfm = {@(x) 1,  @(x) -4*x^3;
             @(x) -1, @(x) 4*x^3};
    case PLANCK_OFF
      fm =  {@(x) x,  @(x) 0;
             @(x) -x, @(x) 0};
      dfm = {@(x) 1,  @(x) 0;
             @(x) -1, @(x) 0};
  endswitch

  data.f = data.df = cell(N, M, N);
  coeff = [n .^ 2 .* mu_a; sigma];
  for i = 1 : N
    for j = 1 : M
      for k = 1 : N
        if (i == 2 && j == M && k == i)
          temp_coeff = omega * rho_b * c_b;
        else
          temp_coeff = 0;
        endif
        data.f{i, j, k} = @(x) coeff(i, j) * fm{i, k}(x) + temp_coeff * x;
        data.df{i, j, k} = @(x) coeff(i, j) * dfm{i, k}(x) + temp_coeff;
      endfor
    endfor
  endfor

  data.g = zeros(N, grid_info.nodes);
  for nn = 0 : K(M) + 1
    ind = gindex(grid_info, M, nn);
    data.g(2, ind) = omega * rho_b * c_b * theta_b;
  endfor

  guess = [ zeros(1, grid_info.nodes) ; ...
    theta_ext + (theta_b - theta_ext) * (0 : (grid_info.nodes - 1)) / (grid_info.nodes - 1) ];
  tol = 1e-6;
  switch (planck)
    case PLANCK_ON
      max_iter = Inf;
    case PLANCK_OFF
      max_iter = 1;
  endswitch

  sol = solve_bvp1d(grid_info, data, guess, tol, max_iter);
  phi = sol(1, :);
  theta = sol(2, :);

  theta_val = get_func_layers(grid_info, theta);
  phi_val = get_func_layers(grid_info, phi);

  temp = intens = Q = cell(1, M);
  # calculate the absorbed energy density
  switch (planck)
    case PLANCK_ON
      for j = 1 : M    
        Q{j} = sigma(j) * (phi_val{j} .- theta_val{j} .^ 4);
      endfor
    case PLANCK_OFF
      for j = 1 : M    
        Q{j} = sigma(j) * phi_val{j};
      endfor
  endswitch
  # calculate the denormalized temperature and radiative intensity
  for j = 1 : M
    temp{j} = theta_val{j} - Kelvin0;
    intens{j} = phi_val{j} / I_norm_coeff * n(j) ^ 2;
  endfor    
endfunction