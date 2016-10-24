# A complex heat transfer problem for 3 layers with various refractive indices

clear all
more off
format long

input_data

M = 3;
n1_range = 1 : 0.1 : 2;
cnt = 0;
for n1_value = n1_range
  printf("n1 = n3 = %f\n", n1_value)
  layer_data(1).n = layer_data(3).n = n1_value;
  [grid_info, theta, phi] = trilayer_calc(layer_data, boundary_data, true);
  [grid_info, theta_wo, phi_wo] = trilayer_calc(layer_data, boundary_data, false);
  xgrid = get_grid_layers(grid_info);
  xgrid_all = concat_cells(xgrid, 1, M);
  theta_all = concat_cells(theta, 1, M);
  theta_wo_all = concat_cells(theta_wo, 1, M);
  rms(++cnt) = sqrt(trapz(xgrid_all, (theta_all .- theta_wo_all) .^ 2));
endfor

frms = fopen("rms/rms.txt", "wt");
for i = 1 : cnt
  fprintf(frms, "%.12f %.12f\n", n1_range(i), rms(i));
endfor
fclose(frms);

figure
plot(n1_range, rms);
xlabel("Relative refractive index");
ylabel("Root mean square deviation");
saveas(1, "rms/rms.eps");
