X_ = csvread("./data/x.csv");
X_ = sort(X_);
n_max = length(X_);

printf("Введите объем выборки (от 1 до %d включительно):\n\t", n_max);
n = str2double(input("N = ", "s"));

X = X_(1:n);
mu = sum(X) / n;
s_square = sum((X - mu).^2) / (n - 1);

printf("Введите уровень доверия (от 0 до 1 включительно):\n\t");
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

N = 1:n;

mu_N = zeros(n, 1);
mu_low_N = zeros(n, 1);
mu_high_N = zeros(n, 1);

for n_i = 1:n
    t_i = tinv((1 + gamma) / 2, n_i - 1);
    mu_N(n_i) = sum(X(1:n_i)) / n_i;
    mu_low_N(n_i) = mu_N(n_i) - s * t_i / sqrt(n_i);
    mu_high_N(n_i) = mu_N(n_i) + s * t_i / sqrt(n_i);
endfor;

figure;
p1 = plot(N, mu, "Color", "black", "LineStyle", "-", "linewidth", 2);
hold on;
p2 = plot(N, mu_N, "Color", "green", "LineStyle", "-.", "linewidth", 2);
p3 = plot(N, mu_low_N, "Color", "red", "LineStyle", "--", "linewidth", 2);
p4 = plot(N, mu_high_N, "Color", "blue", "LineStyle", ":", "linewidth", 2);
##h = plot(N,mu, N,mu_N, N,mu_low_N, N,mu_high_N);
hold off;
##legend(h, {"mu N", "mu n", "mu low n", "mu high n"}, "location", "southeast");
xlabel("Объем выборки n");
grid on;

s_square_N = zeros(n, 1);
sigma_low_N = zeros(n, 1);
sigma_high_N = zeros(n, 1);

for n_i = 1:n
    s_square_N(n_i) = sum((X(1:n_i) - mu).^2) / (n_i - 1);
    h_low_i = chi2inv((1 + gamma) / 2, n_i - 1);
    h_high_i = chi2inv((1 - gamma) / 2, n_i - 1);
    sigma_low_N(n_i) = (n_i - 1) * s_square_N(n_i) / h_low_i;
    sigma_high_N(n_i) = (n_i - 1) * s_square_N(n_i) / h_high_i;
endfor;

figure;
p1 = plot(N, s_square, "Color", "black", "LineStyle", "-", "linewidth", 2);
hold on;
p2 = plot(N, s_square_N, "Color", "green", "LineStyle", "-.", "linewidth", 2);
p3 = plot(N, sigma_low_N, "Color", "red", "LineStyle", "--", "linewidth", 2);
p4 = plot(N, sigma_high_N, "Color", "blue", "LineStyle", ":", "linewidth", 2);
##h = plot(N,mu, N,mu_N, N,mu_low_N, N,mu_high_N);
hold off;
##legend(h, {"mu N", "mu n", "mu low n", "mu high n"}, "location", "southeast");
xlabel("Объем выборки n");
grid on;
