warning("off", "Octave:data-file-in-path");
X = csvread('x.csv');

X = sort(X);
n = length(X);

M_max = max(X);
M_min = min(X);

fprintf("а) вычисление максимального значения M_max и\
		минимального значения M_min:\n");
fprintf("\tM_max = %+.4f\n\tM_min = %+.4f\n", M_max, M_min);

R = M_max - M_min;

fprintf("\nб) вычисление размаха R выборки:\n");
fprintf("\tR = %.4f\n", R);

Mu = sum(X) / n;
S_quad = sum((X - Mu) .^2) / (n - 1);
S = sqrt(S_quad);

fprintf("\nв) вычисление оценок Mu и S_quad математического\
		ожидания MX и дисперсии DX:\n");
fprintf("\tMu = %.4f\n\tS_quad = %.4f\n", Mu, S_quad);

fprintf("\nг) группировка значений выборки в m = [log2(n)] + 2\
		интервала:\n");

m = floor(log2(n)) + 2;
fprintf("\tКоличество интервалов m = %d\n\n", m);
delta = (X(n) - X(1)) / m;
J_limits = M_min : delta : M_max;
ni = zeros(m, 1);
