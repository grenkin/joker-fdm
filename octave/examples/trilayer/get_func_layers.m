function func_layers = get_func_layers (grid_info, func)
  func_layers = cell(1, grid_info.M);
  Lcum = 0;
  for j = 1 : grid_info.M
    func_layers{j} = func(grid_info.first_index(j) + (0 : grid_info.K(j)));
    Lcum += grid_info.L(j);
  endfor
endfunction
