X_ = csvread("./data/x.csv");
X_ = sort(X_);
n_max = length(X_);

n = n_max;

begin = 1;

X = X_(1:n);
mu = sum(X) / n;

if (n > 1)
    s_square = sum((X - mu).^2) / (n - 1);
else
    s_square = 0;
endif

printf("Введите уровень доверия (вещественное от 0 до 1 включительно):\n\t");
gamma = str2double(input("gamma = ", "s"));

printf("\nТочечные оценки MX и DX:\n");
printf("\tmu       = %+.4f\n\tS_square = %+.4f\n", mu, s_square);

s = sqrt(s_square);

t = tinv((1 + gamma) / 2, n - 1);
n_sqrt = sqrt(n);

mu_low = mu - s * t / n_sqrt;
mu_high = mu + s * t / n_sqrt;

printf("\nНижняя и верхняя границы для gamma-доверительного интервала для\
        \nматематического ожидания MX:\n");
printf("\tmu_low  = %+.4f\n\tmu_high = %+.4f\n", mu_low, mu_high);

h_low = chi2inv((1 + gamma) / 2, n - 1);
h_high = chi2inv((1 - gamma) / 2, n - 1);
sigma_low = (n - 1) * s_square / h_low;
sigma_high = (n - 1) * s_square / h_high;

printf("\nНижняя и верхняя границы для gamma-доверительного интервала для\
        \nматематического ожидания DX:\n");
printf("\tsigma_low  = %+.4f\n\tsigma_high = %+.4f\n", sigma_low, sigma_high);

mu_n = zeros(n - begin + 1, 1) + mu;
mu_N = zeros(n - begin + 1, 1);
mu_low_N = zeros(n - begin + 1, 1);
mu_high_N = zeros(n - begin + 1, 1);

for n_i = begin:n
    t_i = tinv((1 + gamma) / 2, n_i - 1);
    mu_N(n_i - begin + 1) = sum(X(1:n_i)) / n_i;
    mu_low_N(n_i - begin + 1) = mu_N(n_i - begin + 1) - s * t_i / sqrt(n_i);
    mu_high_N(n_i - begin + 1) = mu_N(n_i - begin + 1) + s * t_i / sqrt(n_i);
endfor;

figure;
p1 = plot((begin:n), mu_n, "Color", "black", "linestyle", "-", "linewidth", 2);
hold on;
p2 = plot((begin:n), mu_N, "Color", "magenta", "linestyle", "--", "linewidth", 2);
p3 = plot((begin:n), mu_low_N, "Color", "red", "linestyle", "-.", "linewidth", 2);
p4 = plot((begin:n), mu_high_N, "Color", "blue", "linestyle", ":", "linewidth", 2);
hold off;
xlabel("Объем выборки n");
grid on;

s_square_n = zeros(n - begin + 1, 1) + s_square;
s_square_N = zeros(n - begin + 1, 1);
sigma_low_N = zeros(n - begin + 1, 1);
sigma_high_N = zeros(n - begin + 1, 1);

for n_i = begin:n
    s_square_N(n_i - begin + 1) = sum((X(1:n_i) - mu).^2) / (n_i - 1);
    h_low_i = chi2inv((1 + gamma) / 2, n_i - 1);
    h_high_i = chi2inv((1 - gamma) / 2, n_i - 1);
    sigma_low_N(n_i - begin + 1) = (n_i - 1) * s_square_N(n_i - begin + 1) / h_low_i;
    sigma_high_N(n_i - begin + 1) = (n_i - 1) * s_square_N(n_i - begin + 1) / h_high_i;
endfor;

figure;
p1 = plot((begin:n), s_square_n, "Color", "black", "linestyle", "-", "linewidth", 2);
hold on;
p2 = plot((begin:n), s_square_N, "Color", "magenta", "linestyle", "--", "linewidth", 2);
p3 = plot((begin:n), sigma_low_N, "Color", "red", "linestyle", "-.", "linewidth", 2);
p4 = plot((begin:n), sigma_high_N, "Color", "blue", "linestyle", ":", "linewidth", 2);
hold off;
xlabel("Объем выборки n");
grid on;


