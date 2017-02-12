function xgrid = get_grid_layers (grid_info)
  xgrid = cell(1, grid_info.M);
  Lcum = 0;
  for j = 1 : grid_info.M
    xgrid{j} = Lcum + (0 : grid_info.K(j)) * grid_info.h(j);
    Lcum += grid_info.L(j);
  endfor
endfunction