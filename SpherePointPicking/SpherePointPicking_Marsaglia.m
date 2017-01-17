%Sphere Point Picking:Marsaglia
x1 = random('uniform',-1,1,1,100);
x2 = random('uniform',-1,1,1,100);
x = zeros(3,1);
y = zeros(3,1);
z = zeros(3,1);
for i=1:1:100
    if x1(i)^2 + x2(i)^2 <= 1
        x(i) = 2*x1(i)*sqrt(1-x1(i)^2-x2(i)^2);
        y(i) = 2*x2(i)*sqrt(1-x1(i)^2-x2(i)^2);
        z(i) = 1-2*(x1(i)^2+x2(i)^2);
    end
end
plot3(x,y,z,'.r');