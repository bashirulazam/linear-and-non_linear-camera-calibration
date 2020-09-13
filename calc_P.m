function P = calc_P(Y,Z)

p1 = Y(1:3,1);
p14 = Y(4,1);
p2 = Y(5:7,1);
p24 = Y(8,1);
p34 = Y(9,1);
p3 = Z;
P = [ p1' p14; p2' p24; p3' p34];


