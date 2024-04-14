step = S / 1000;
x_coords = (M_min - R) : step : (M_max + R);
f_density_normal = normpdf(x_coords, Mu, S);

figure
hold on;
bar(J_middles, fn_values, 1, 'y')
plot(x_coords, f_density_normal, 'b', 'LineWidth', 2);
grid on;
hold off;

fprintf("\nе) построение на другой координатной плоскости\
		графика \n   эмпирической функции распределения и\
		функции \n   распределения нормальной случайной ве\
		личины с \n   математическим ожиданием Mu и диспер\
		сией S_quad\n"); 

t_arr = zeros(n + 2, 1);
t_arr(1) = X(1) - 1;
t_arr(n + 2) = X(n) + 1;

for i = 2 : n + 1
    t_arr(i) = X(i - 1);
endfor

f_emp = zeros(length(t_arr), 1);

for i = 1 : length(t_arr)
    count = 0;
    
    for j = 1 : n
        if (t_arr(i) >= X(j))
            count++;
        end
    endfor
    
    f_emp(i) = count / n;
endfor
