# Calculate the coefficient G_i
# n1 -- refractive index in the domain to the left (n_i)
# n2 -- refractive index in the domain to the right (n_{i+1})
# Returned value -- G_i
# n1, n2 are positive or equal Inf
function retval = calc_G (n1, n2)
  if (n1 == n2)
    retval = Inf;
    return;
  endif
  int1 = quadcc(integrand_G_numer(n1, n2), 0, 1);
  int2 = quadcc(integrand_G_denom(n1, n2), 0, 1);
  retval = (n1 ^ 2 * int1) / (3 * int2);
endfunction

function retfunc = integrand_G_numer (n1, n2)
  n12 = n1 / n2;
  retfunc = @(mu_vect) arrayfun(@(mu)(mu * (1 - R_ij(mu, n12))), mu_vect);
endfunction

function retfunc = integrand_G_denom (n1, n2)
  n12 = n1 / n2;
  n21 = n2 / n1;
  retfunc = @(mu_vect) arrayfun(@(mu)(mu ^ 2 * (R_ij(mu, n12) + R_ij(mu, n21))), mu_vect);
endfunction
