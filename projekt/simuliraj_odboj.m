function [pos_new, v_new, index_novega_clena, xs, ys] = simuliraj_odboj(pos, v, veriznica, g, koeficient_izgube, casovni_interval)
%function [pos_new, v_new, index_novega_clena] = simuliraj_odboj(pos, v, veriznica, g, koeficient_izgube, casovni_interval)
%   Detailed explanation goes here
x = pos(1,1);
y = pos(2,1);
vx = v(1,1);
vy = v(2,1);
N = size(veriznica, 2) - 1;

%plot([pos(1,1), pos(1,1)+v(1,1)], [pos(2,1), pos(2,1)+v(2,1)])
% enačbe premic y = kx + n
y2 = veriznica(2, 2:end);
y1 = veriznica(2, 1:end-1);
x2 = veriznica(1, 2:end);
x1 = veriznica(1, 1:end-1);
k = (y2 - y1) ./ (x2 - x1);
n = y1 - k .* x1;

if abs(vx) < 1e-6
    pos_new = pos;
    seka_palico = (x1 <= x) & (x <= x2);
    assert(sum(seka_palico) >= 1, 'kroglica gre ven iz verižnice');
    [val, index_novega_clena] = max(seka_palico);
    koef = k(1, index_novega_clena);
    pos_new(2,1) = y1(1, index_novega_clena) + koef .* (x - x1(1,index_novega_clena));
    delta_h = -pos(2,1) + pos_new(2,1);
    a = g./2;
    b = -vy;
    c = delta_h;
    intersection_time = (-b + sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
else
   
    % parametri parabole leta
    a = -g ./ (2.*vx.^2);
    b = vy./vx + (x.*g)./(vx.^2);
    c = y - (g.*x.^2)./(2.*vx.^2) - (vy.*x) ./ vx;
    
    yf = @(x) a.*x.^2 + b .*x + c;
    
    %xs = linspace(0, 30, 1000);
    %ys = yf(xs);
    %plot(xs, ys, 'g')
    
    % presečišča parabole in premic
    a = a .* ones(1, N);
    b = b - k;
    c = c - n;
    
    D = sqrt(b.^2 - 4.*a.*c);
    % prva rešitev <= druga rešitev, ker a negativen
    x_solution1 = (-b + D) ./ (2.*a);
    x_solution2 = (-b - D) ./ (2.*a);
    time1 = (x_solution1 - x) ./ vx;
    time2 = (x_solution2 - x) ./ vx;
    
    % če sta oba pozitivna, potem vzamemo manjši čas
    time1(time1<1e-5) = inf;
    time2(time2<1e-5) = inf;
    time = min(time1, time2);
    x_solution = time .* vx + x;
    
    % pogledamo, katere rešitve so med spojema palice
    seka_palico = (x1 <= x_solution) & (x_solution <= x2);
    assert (sum(seka_palico) >= 1, "kroglica gre ven iz verižnice");
    
    % med tistimi, ki sekajo palico, izberemo tisto, ki ima najmanjši čas
    relevant_indices = find(seka_palico == 1);
    relevant_values = time(relevant_indices);
    [intersection_time, min_index_local] = min(relevant_values);
    min_index_original = relevant_indices(min_index_local);
    index_novega_clena = min_index_original;

    pos_new = [x_solution(1, index_novega_clena);
            yf(x_solution(1, index_novega_clena))];
end

% hitrost pred trkom
v = [vx;
     vy - g.* intersection_time];

k = k(1, index_novega_clena); % gledamo za tisti člen, od katerega se bo odbilo
n = [-k; 1];               % normalni vektor
n = n / norm(n);           % enotski normalni vektor
t = [1; k];
t = t / norm(t);           % enotski tangencialni vektor
vn = dot(v, n) * n;        % projekcija na normalo
vt = dot(v, t) * t;        % projekcija na tangento

v_new = vt - vn;           % nova hitrost po odboju
v_new = v_new .* koeficient_izgube;

xf = @(t) x + vx .* t;
yf = @(t) y + vy .* t - g .* t.^2 ./ 2;
t = 0 : casovni_interval : intersection_time;
xs = xf(t);
ys = yf(t);

end

