N = 21;
n = (N-1)/2;
a = 0.25;
B = 2;

c = (B - a) ./ (2.*a);

function value = f_helper(theta1, c, n)
    tan_theta1 = tan(theta1);
    tan_thetas = (1:n) .* tan_theta1;
    cos_thetas = 1 ./ (sqrt(1 + tan_thetas.^2))

    value = sum(cos_thetas) - c;
end

f = @(theta1) f_helper(theta1, c, n);
theta1 = fsolve(f, 0.5)

% Compute angles
tan_theta1 = tan(theta1);
tan_thetas = (1:n) .* tan_theta1;
cos_thetas = 1 ./ (sqrt(1 + tan_thetas.^2))
sin_thetas = tan_thetas .* cos_thetas;


x0 = a/2;
y0 = 0;
X = [x0];
Y = [y0];
x = x0;
y = y0;
for i=1:n
    x = x + a .* cos_thetas(1, i);
    y = y + a .* sin_thetas(1, i);
    X = [X, x];
    Y = [Y, y];
end

X = [-flip(X), X];
Y = [flip(Y), Y];

plot(X, Y, '-o')
axis equal
hold on

obesisceL = [X(1,1);Y(1,1)];
obesisceD = [X(1,end);Y(1,end)];
L = ones(1, N) .* a;
M = L;
x = diskrVeriznica([-1; -1], obesisceL, obesisceD, L, M);

Y - x(2, :)