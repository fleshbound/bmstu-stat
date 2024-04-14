xs = (M_min - R) : step : (M_max + R);
f_norm = normcdf(x_coords, Mu, S);

figure
hold on;
plot(xs, f_norm, ':r', 'linewidth', 2);
stairs(t_arr, f_emp, 'b', 'linewidth', 2);
grid on;
hold off;
