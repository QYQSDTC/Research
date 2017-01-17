%Sphere Point Picking:Normal
r = 1;
u = random('uniform',0,1,1,100);
v = random('uniform',0,1,1,100);
theta = 2*pi*u;
phi = acos(2*v-1);
x = zeros(3,1);
y = zeros(3,1);
z = zeros(3,1);
for i=1:1:100
x(i) = r*sin(phi(i))*cos(theta(i));
y(i) = r*sin(phi(i))*sin(theta(i));
z(i) = r*cos(phi(i));
end
plot3(x,y,z,'.r');