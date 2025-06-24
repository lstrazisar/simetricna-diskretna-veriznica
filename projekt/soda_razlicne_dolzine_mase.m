N = 10;      % skupno število členov
n = (N)/2;  % število členov na eni polovici
L = [3, 7, 2, 4, 4]; % dolžine členov
M = [2, 5, 1, 4, 5];      % mase členov
B = 30; % razdalja med krajiščema

c = B/2; % polovica razdalje med krajiščema

function value = f_helper(theta1, c, n, L, M)
    tan_theta1 = tan(theta1);
    tan_thetas = zeros(1, n);
    for k=1:n
        koeficient = (2 .* (sum(M(1, 1:k)) .* L(1,k) - M(1,k).*L(1,k).*0.5)) ./ (M(1,1).*L(1,k));
        tanges_k = koeficient .* tan_theta1;
        tan_thetas(1,k) = tanges_k;
    end
    cos_thetas = 1 ./ (sqrt(1 + tan_thetas.^2));
    vsota = sum(cos_thetas .* L)
    value = sum(cos_thetas .* L) - c;
end

f = @(theta1) f_helper(theta1, c, n, L, M);
[theta1, val] = fsolve(f, 0.5)
theta1 = abs(theta1);

% Compute angles
tan_theta1 = tan(theta1);
tan_thetas = zeros(1, n);
for k=1:n
    koeficient = (2 .* (sum(M(1, 1:k)) .* L(1,k) - M(1,k).*L(1,k).*0.5)) ./ (M(1,1).*L(1,k));
    tanges_k = koeficient .* tan_theta1;
    tan_thetas(1,k) = tanges_k;
end
cos_thetas = 1 ./ (sqrt(1 + tan_thetas.^2));
sin_thetas = tan_thetas .* cos_thetas;


x0 = 0;
y0 = 0;
X = [x0];
Y = [y0];
x = x0;
y = y0;
for i=1:n
    x = x + L(1,i) .* cos_thetas(1, i);
    y = y + L(1,i) .* sin_thetas(1, i);
    X = [X, x];
    Y = [Y, y];
end

X = [-flip(X(1, 2:end)), X];
Y = [flip(Y(1, 2:end)), Y];

plot(X, Y, '-o')
axis equal
hold on

L = [flip(L), L];
M = [flip(M), M];
obesisceL = [X(1,1);Y(1,1)];
obesisceD = [X(1,end);Y(1,end)];
resitev = diskrVeriznica([-1; -1], obesisceL, obesisceD, L, M);

Y - resitev(2, :)

function energy = calculate_potential_energy(x, y, L, M)
    x_mid = (x(1, 1:end-1) + x(1, 2:end)) ./ 2;
    y_mid = (y(1, 1:end-1) + y(1, 2:end)) ./ 2;
    energy = sum(M .* y_mid);
end

energy_with_angles = calculate_potential_energy(X ,Y ,L, M)
energy_regular = calculate_potential_energy(resitev(1,:), resitev(2,:), L, M)
razlika = energy_with_angles - energy_regular