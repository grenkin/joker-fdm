# Calculate the coefficient G
# n1, n2 > 0
function retval = calc_G (n1, n2)
  if (n1 == n2)
    retval = Inf;
    return;
  endif
  int1 = quadcc(integrand1(n1, n2), 0, 1);
  int2 = quadcc(integrand2(n1, n2), 0, 1);
  retval = (n1 ^ 2 * int1) / (3 * int2);
endfunction

function retfunc = integrand1 (n1, n2)
  n12 = n1 / n2;
  retfunc = @(mu_vect) arrayfun(@(mu)(mu * (1 - R_ij(mu, n12))), mu_vect);
endfunction

function retfunc = integrand2 (n1, n2)
  n12 = n1 / n2;
  n21 = n2 / n1;
  retfunc = @(mu_vect) arrayfun(@(mu)(mu ^ 2 * (R_ij(mu, n12) + R_ij(mu, n21))), mu_vect);
endfunction

function retval = R_ij (mu, nij)
  if (mu == 0)
    retval = 1;
  else
    frac1 = (psi_ij(mu, nij) - nij * mu) / (psi_ij(mu, nij) + nij * mu);
    frac2 = (nij * psi_ij(mu, nij) - mu) / (nij * psi_ij(mu, nij) + mu);
    retval = 0.5 * (frac1 ^ 2 + frac2 ^ 2);
  endif
endfunction

function retval = psi_ij (mu, nij)
  radicand = 1 - nij ^ 2 * (1 - mu ^ 2);
  if (radicand > 0)
    retval = sqrt(radicand);
  else
    retval = 0;
  endif
endfunction
