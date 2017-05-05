1;

function retval = gamma_func (emiss)
  retval = emiss / ( 2 * (2 - emiss) );
endfunction

function retval = calc_layer (l)
  alpha = 1 / (3 - l.A * l.omega);
  retval.B = l.n ^ 2 * alpha;
  retval.C = l.Nc * l.n ^ 2;
  retval.K = l.n ^ 2 * (1 - l.omega);
  retval.L = l.L;
endfunction

function retval = calc_boundary (b)
  emiss = 1 - b.R;
  gam = gamma_func(emiss);
  retval.tilde_gamma = b.n ^ 2 * gam;
  retval.thetab = b.thetab;
endfunction
