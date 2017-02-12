function ret = print_func_layers (filename, xgrid_layers, func_layers)
  ret = 0;
  f = fopen(filename, "wt");
  M = length(xgrid_layers);
  for j = 1 : M
    for i = 1 : length(xgrid_layers{j})
      fprintf(f, "%.12f   %.12f\n", xgrid_layers{j}(i), func_layers{j}(i));
    endfor
    fprintf(f, "\n\n");
  endfor
  fclose(f);
endfunction