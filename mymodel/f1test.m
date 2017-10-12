
g_a = 0.03;
g_f = 20;

R = 1/0.8;

L = 2;
p = 50;



%L = R*2*start_corner;

delta = L/p;

a1 = 0:delta:L;
a2 = 0.025;

ar = (pi*R+L)/(2*R) - a1/R;

x=(R + a2 + g_a.*cos(g_f.*ar)).*cos(ar);
y=(R + a2 + g_a.*cos(g_f.*ar)).*sin(ar);

%ar_grad = ar*180/pi
figure
plot(x,y, 'g');
hold on;
alpha1=a1(1:4:end);
alpha2=a2(1:4:end);

q=1+alpha2./R;
ar = (pi*R+L)/(2*R) - alpha1./R;
w = q+g_a*K*cos(g_f*ar);
z = g_a*g_f*K*sin(g_f*ar);


x=(R + a2 + g_a.*cos(g_f.*ar)).*cos(ar)
y=(R + a2 + g_a.*cos(g_f.*ar)).*sin(ar)

r11=(w.*sin(ar)+z.*cos(ar))
r12=(-w.*cos(ar)+z.*sin(ar))

r21=cos(ar)
r22=sin(ar)

quiver (x,y, r11, r12, 'r')
hold on;
quiver (x,y, r21,r22,'b')


%============================================
start_corner = L/(2*R);
delta = 2 * start_corner / p
psi = -start_corner:delta:start_corner;


x=(R + g_a.*cos(g_f.*psi)).*cos(psi);
y=(R + g_a.*cos(g_f.*psi)).*sin(psi);

%plot(x,y);


 %clf;
 [x,y] = meshgrid (1:2:20);
x = -0.87609
y =  0.85087
u =  0.77983
v =  0.60459
 
 
 %h = quiver (x,y, u, v);
