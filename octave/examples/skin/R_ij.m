# Calculate the reflection coefficient by Fresnel formula
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
