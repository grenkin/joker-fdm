# Weighted sums for different wavelengths
# in the absence and in the presence of nanoparticles
# with taking into account refraction and reflection at the interfaces
# without the term theta^4

clear all
more off

calc

freqs = [FREQ_UV, FREQ_VIS, FREQ_IR];
freqs_num = length(freqs);
global freqs_w = [0.05, 0.5, 0.45];

temp = intens = Q = cell(1, freqs_num);
temp_nano = intens_nano = Q_nano = cell(1, freqs_num);
for fi = 1 : freqs_num
  freq = freqs(fi);
  [grid_info3, temp{fi}, intens{fi}, Q{fi}] = skin_calc(freq, PLANCK_OFF, REFR_ON, NANO_OFF);
  [grid_info4, temp_nano{fi}, intens_nano{fi}, Q_nano{fi}] = skin_calc(freq, PLANCK_OFF, REFR_ON, NANO_ON);
endfor

function ret = get_sup (v)
  global freqs_w
  M = length(v{1});
  ret = cell(1, M);
  for j = 1 : M
    ret{j} = v{1}{j} * freqs_w(1) + v{2}{j} * freqs_w(2) + v{3}{j} * freqs_w(3);  
  endfor
endfunction

temp_sup = get_sup(temp);
intens_sup = get_sup(intens);
Q_sup = get_sup(Q);
temp_nano_sup = get_sup(temp_nano);
intens_nano_sup = get_sup(intens_nano);
Q_nano_sup = get_sup(Q_nano);

xgrid3 = get_grid_layers(grid_info3);
xgrid4 = get_grid_layers(grid_info4);
L_total = sum(grid_info3.L);
xgrid3_all = concat_cells(xgrid3, 1, 3);
xgrid4_all = concat_cells(xgrid4, 1, 4);



figure(1)
plot(1e6 * xgrid3_all, concat_cells(temp_sup, 1, 3), "b", ...
  1e6 * xgrid4_all, concat_cells(temp_nano_sup, 1, 4), "m");
xlabel('z, micrometers')
ylabel('Temperature, C')
#xlim([0, 1e6 * L_total])
# ylim([36, 38.5])
filename = "4/temp_sup";
saveas(1, [filename ".jpg"]);
print_func_layers([filename ".txt"], xgrid3, temp_sup);
print_func_layers([filename "_nano.txt"], xgrid4, temp_nano_sup);

figure(2)
plot(1e6 * concat_cells(xgrid3, 1, 2), 1e-6 * concat_cells(Q_sup, 1, 2), "b", ...
  1e6 * concat_cells(xgrid4, 1, 3), 1e-6 * concat_cells(Q_nano_sup, 1, 3), "m");
xlabel("z, micrometers");
ylabel("Absorbed energy density Q, W / cm^3");
filename = "4/Q_sup";
saveas(2, [filename ".jpg"]);
print_func_layers([filename ".txt"], xgrid3, Q_sup);
print_func_layers([filename "_nano.txt"], xgrid4, Q_nano_sup);

figure(3)
plot(1e6 * concat_cells(xgrid3, 1, 2), concat_cells(intens_sup, 1, 2), "b", ...
  1e6 * concat_cells(xgrid4, 1, 3), concat_cells(intens_nano_sup, 1, 3), "m");
xlabel("z, micrometers");
ylabel('Radiative intensity I_0, W / (m^2 * srad)', "interpreter", "tex");
filename = "4/intens_sup";
saveas(3, [filename ".jpg"]);
print_func_layers([filename ".txt"], xgrid3, intens_sup);
print_func_layers([filename "_nano.txt"], xgrid4, intens_nano_sup);


# graphs for all frequences

figure(4)
plot(1e6 * concat_cells(xgrid3, 1, 2), 1e-6 * freqs_w(1) * concat_cells(Q{1}, 1, 2), "b", ...
  1e6 * concat_cells(xgrid3, 1, 2), 1e-6 * freqs_w(2) * concat_cells(Q{2}, 1, 2), "g", ...
  1e6 * concat_cells(xgrid3, 1, 2), 1e-6 * freqs_w(3) * concat_cells(Q{3}, 1, 2), "r", ...
  1e6 * concat_cells(xgrid3, 1, 2), 1e-6 * concat_cells(Q_sup, 1, 2), "k");
xlabel("z, micrometers");
ylabel("Absorbed energy density Q, W / cm^3");
saveas(4, "4/Q_all_freqs.jpg");

figure(5)
plot(1e6 * concat_cells(xgrid4, 1, 2), 1e-6 * freqs_w(1) * concat_cells(Q_nano{1}, 1, 2), "b", ...
  1e6 * concat_cells(xgrid4, 1, 2), 1e-6 * freqs_w(2) * concat_cells(Q_nano{2}, 1, 2), "g", ...
  1e6 * concat_cells(xgrid4, 1, 2), 1e-6 * freqs_w(3) * concat_cells(Q_nano{3}, 1, 2), "r", ...
  1e6 * concat_cells(xgrid4, 1, 2), 1e-6 * concat_cells(Q_nano_sup, 1, 2), "k");
xlabel("z, micrometers");
ylabel("Absorbed energy density Q, W / cm^3");
saveas(5, "4/Q_all_freqs_nano.jpg");
