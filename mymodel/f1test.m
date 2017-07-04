
g_a = 0.03;
g_f = 20;

R = 1/0.8;

L = 2;
p = 50;



%L = R*2*start_corner;

delta = L/p;

a1 = 0:delta:L;
a2 = 0;

ar = (pi*R+L)/(2*R) - a1/R;

x=(R + a2 + g_a.*cos(g_f.*ar)).*cos(ar);
y=(R + a2 + g_a.*cos(g_f.*ar)).*sin(ar);

%ar_grad = ar*180/pi

plot(x,y);


%============================================
start_corner = L/(2*R);
delta = 2 * start_corner / p
psi = -start_corner:delta:start_corner;


x=(R + g_a.*cos(g_f.*psi)).*cos(psi);
y=(R + g_a.*cos(g_f.*psi)).*sin(psi);

%plot(x,y);



