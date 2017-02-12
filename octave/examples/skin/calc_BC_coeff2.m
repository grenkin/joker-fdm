# Calculate the coefficient in the boundary condition at the left point
# n1 > 0 -- refractive index in the first domain
function retval = calc_BC_coeff2 (n1)
  int1 = quadcc(integrand_BC_numer(n1), 0, 1);
  int2 = quadcc(integrand_BC_denom(n1), 0, 1);
  retval = (1 / 3 - int1) / (1 + 3 * int2);
endfunction

function retfunc = integrand_BC_numer (n1)
  n0 = 1;
  n01 = n0 / n1;
  retfunc = @(mu_vect) arrayfun(@(mu)(mu ^ 2 * R_ij(mu, n01)), mu_vect);
endfunction

function retfunc = integrand_BC_denom (n1)
  n0 = 1;
  n10 = n1 / n0;
  retfunc = @(mu_vect) arrayfun(@(mu)(mu ^ 2 * R_ij(mu, n10)), mu_vect);
endfunction
