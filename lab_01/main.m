warning("off", "Octave:data-file-in-path");
X = csvread('x.csv');
X = sort(X);
n = length(X);

fprintf("а) вычисление максимального значения M_max и минимального значения \
M_min:\n");
M_max = max(X);
M_min = min(X);
fprintf("\tM_max = %+.4f\n\tM_min = %+.4f\n", M_max, M_min);

fprintf("\nб) вычисление размаха R выборки:\n");
R = M_max - M_min;
fprintf("\tR = %.4f\n", R);

fprintf("\nв) вычисление оценок Mu и S_square математического ожидания MX и \
дисперсии DX:\n");
Mu = sum(X) / n;
S_square = sum((X - Mu) .^2) / (n - 1);
fprintf("\tMu = %.4f\n\tS_square = %.4f\n", Mu, S_square);

S = sqrt(S_square);

fprintf("\nг) группировка значений выборки в m = [log2(n)] + 2 интервала:\n");
m = floor(log2(n)) + 2;
fprintf("\tКоличество интервалов m = %d\n\n", m);
delta = (X(n) - X(1)) / m;
J_limits = M_min : delta : M_max;
ni = zeros(m, 1);

for i = 1 : m
    count = 0;
    
    for x = X
        if (i == m) && (x >= J_limits(i)) && (x <= J_limits(i + 1))
            count++;
        elseif (x >= J_limits(i)) && (x < J_limits(i + 1))
            count++;
        endif
    endfor
    
    if (i == m)
        fprintf("\t%d) [%+.3f; %+.3f), n%d = %d\n", i, J_limits(i), 
                J_limits(i + 1), i, count);
    else
        fprintf("\t%d) [%+.3f; %+.3f], n%d = %d\n", i, J_limits(i),
                J_limits(i + 1), i, count);
    endif
    
    ni(i) = count;
endfor

fprintf("\nд) построение на одной координатной плоскости гистограммы\n\
   и графика функции плотности распределния вероятностей\n   \
   нормальной случайной величины с математическим\n   \
   ожиданием Mu и дисперсией S_square\n");

J_middles = zeros(m, 1);

for i = 1 : m
    J_middles(i) = (J_limits(i) + J_limits(i + 1)) / 2;
endfor

fn_values = zeros(m, 1);

for i = 1 : m
    fn_values(i) = ni(i) / (n * delta);
endfor

step = 0.001;
x_coords = (M_min - 1) : step : (M_max + 1);
f_density_normal = normpdf(x_coords, Mu, S);

figure
hold on;
bar(J_middles, fn_values, 1, 'y')
plot(x_coords, f_density_normal, 'b', 'LineWidth', 2);
grid on;
hold off;

fprintf("\nе) построение на другой координатной плоскости графика \n\
    эмпирической функции распределения и функции\n\
    распределения нормальной случайной величины с\n\
    математическим ожиданием Mu и дисперсией S_square\n"); 

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

xs = (M_min - 1) : step : (M_max + 1);
f_norm = normcdf(x_coords, Mu, S);

figure
hold on;
plot(xs, f_norm, ':r', 'linewidth', 2);
stairs(t_arr, f_emp, 'b', 'linewidth', 2);
grid on;
hold off;