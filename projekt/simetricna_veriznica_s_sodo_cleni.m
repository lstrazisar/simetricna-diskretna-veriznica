function veriznica = simetricna_veriznica_s_sodo_cleni(obesisceL, obesisceD, L, M, risi)
% function veriznica = simetricna_veriznica_s_sodo_cleni(obesisceL, obesisceD, L, M, risi)
% simetricna_veriznica_s_sodo_cleni reši problem simetrične verižnice s
% sodo mnogo členi: preko f solve najde kot theta1, ki je kot med vodoravno
% premico in najnižjo palico. Nato preko tega kota izračuna še vse ostale
% in vrne koordinate spojev vključno z obesiščema
% N ... število palic
% n := N/2 ... število palic na eni polovici
% vhod:
% obesisceL = [x_0;y_0],
% obesisceD = [x_N;y_N],  obesišči morata biti na enaki višini
% L = dolzine palic levega dela verižnice od leve proti desni (vektor 1xn).
% M = mase palic levega dela verižnice od leve proti desni (vektor 1xn).
% risi = 1, če hočemo izris verižnice, 0 drugače
%
% izhod:
% veriznica je 2x(N+1) matrika, v prvi vrstici so x-koordinate spojev, v 2.
% vrstici so y-koordinate

assert(abs(obesisceL(2,1) - obesisceD(2,1)) < 1e-9, 'Obesišči morata biti na enaki višini')
L = flip(L); % obrnemo, ker bomo delali za desno stran
M = flip(M); % obrnemo, ker bomo delali za desno stran
B = obesisceD(1,1) - obesisceL(1,1); % razdalja med obesiščema v x smeri
c = B/2; % polovica razdalje, ker bomo delali samo za eno stran
n = size(L, 2);
N = 2*n;

assert (sum(L) >= c, "palice so prekratke za verižnico med obesiščema")

f = @(theta1) f_helper(theta1, c, n, L, M);
[theta1, val] = fsolve(f, 0.5);
theta1 = abs(theta1);

% Compute angles
tan_theta1 = tan(theta1);
kumulativne_vsote_mas = cumsum(M);
koeficienti = (2 .* (kumulativne_vsote_mas .* L - M .* L .* 0.5)) ./ (M(1,1).*L);
tan_thetas = koeficienti .* tan_theta1;
cos_thetas = 1 ./ (sqrt(1 + tan_thetas.^2));
sin_thetas = tan_thetas .* cos_thetas;

X = L .* cos_thetas;
Y = L .* sin_thetas;
X = cumsum(X);
Y = cumsum(Y);

% zrcalimo, da dobimo še drugo polovico
X = [-flip(X), 0, X];
Y = [flip(Y), 0, Y];

% premaknemo tako, da je obesišče na pravem mestu
X = X + (obesisceL(1,1)- X(1,1));
Y = Y + (obesisceL(2,1)- Y(1,1));
veriznica = [X;Y];

if risi
    plot(X, Y, '-o', 'Color', 'b')
    axis equal
end

end

% funkcija, katere ničlo iščemo
% izpeljana je iz ravnovesja sil in navorov
function value = f_helper(theta1, c, n, L, M)
    tan_theta1 = tan(theta1);
    kumulativne_vsote_mas = cumsum(M);
    koeficienti = (2 .* (kumulativne_vsote_mas .* L - M .* L .* 0.5)) ./ (M(1,1).*L);
    tan_thetas = koeficienti .* tan_theta1;
    cos_thetas = 1 ./ (sqrt(1 + tan_thetas.^2));
    % hočemo, da je vsota x komponent palic enaka polovici razdalje v x
    % smeri med obesiščema
    value = sum(cos_thetas .* L) - c;
end