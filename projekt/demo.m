obesisceL = [1; 5];
obesisceD = [12; 5];

L = [1, 2, 0.7, 1.4, 2.3];
M = [2, 5, 1, 4, 5]; 

veriznica = simetricna_veriznica_s_sodo_cleni(obesisceL, obesisceD, L, M, 1);
hold on;
axis equal

stevilo_odbojev = 5;
pos = [7.3;6];
v = [5;1];
g = 9.81;
koeficient_izgube = 0.9;
casovni_interval = 0.01;

x = [];
y = [];
for i=1:stevilo_odbojev
    [pos, v, index_novega_clena, xs ,ys] = simuliraj_odboj(pos, v, veriznica, g, koeficient_izgube, casovni_interval);
    x = [x, xs];
    y = [y, ys];
end

% dodamo zadnjo točko
x = [x, pos(1,1)];
y = [y, pos(2,1)];

koncni_clen = index_novega_clena

% vizualizacija
h = animatedline('Color', 'r', 'LineWidth', 1);
for i = 1:size(x,2)
    addpoints(h, x(i), y(i));
    drawnow;
    pause(casovni_interval);  % Adjust speed
end