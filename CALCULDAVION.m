%%%%%%%%%%%%%%%%%%%%%%%%sanssemelle
% e1=1;
% e2=1;
% e3=1;
e1=2.5;
e2=1.5;
e3=2;
b=500;
h=400;
TY=10000;
w=160;
x=60;
S=9600;
Iz1=(b*e3^3)/3+b*e3*h^2+e1*h^3/6+pi*e2*h^3/8;
Iz2=4*(w*(x^3)/12+w*x*(h/2)^2);
Iz=Iz1+Iz2;
%%%%%%%%%%%%%%%%%%%%%%%%
% Y=-200:20:200;
% tao1M = (e1*h^2*(1/6*((2+3*pi*e3/e2*b/h)/(2+pi*e1/e2))+1/8*(1-4*Y.^2/h^2))*TY/Iz1)/e1;
% t1M = (e1*h^2*(1/6*((2+3*pi*(e3*b)/(e2*h)+3*pi*S/(e2*h))/(2+pi*e1/e2))+1/8*(1-4*Y.^2/h^2))*TY/Iz)/e1;
% plot(Y,tao1M,'b.-','LineWidth',2,'MarkerSize',20)
% hold on;
% plot(Y,t1M,'r.-','LineWidth',2,'MarkerSize',20)
% axis([-inf,inf,0,5]);
% set(gca,'XMinorTick','on');
% xlabel('Axis Y (mm)','FontWeight','bold');
% ylabel('La contrainte (MPa)','FontWeight','bold');
% legend('Cas 2(sans semelle)','Cas 4(avec semelle)')
% grid minor;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o=0:1/18*pi:pi;
tao2P = (e2*h^2*(1/3*((3*e3/e1*b/h-1)/(pi+2*e2/e1))+1/4*sin(o))*TY/Iz1)/e2;
t2P = (e2*h^2*(1/3*((3*e3/e1*b/h+3*S/e1/h-1)/(pi+2*e2/e1))+1/4*sin(o))*TY/Iz)/e2;
plot(o,tao2P,'b.-','LineWidth',2,'MarkerSize',20)
hold on;
plot(o,t2P,'r.-','LineWidth',2,'MarkerSize',20)
axis([-inf,inf,0,4]);
% set(gca,'XMinorTick','on');
set(gca, 'xtick', 0:1/9*pi:pi);
set(gca, 'xticklabel', {'0', '¦Ð/9', '2¦Ð/9', '¦Ð/3', '4¦Ð/9', '5¦Ð/9', '2¦Ð/3', '7¦Ð/9', '8¦Ð/9','¦Ð'})
xlabel('LAngle (¡ã)','FontWeight','bold');
ylabel('La contrainte (MPa)','FontWeight','bold');
legend('Cas 2(sans semelle)','Cas 4(avec semelle)')
grid minor;

%% %
% Z=-500:50:500;
% tao3 = (-1/2*e3.*Z*h*TY/Iz1)/e3;
% t3=(-TY/Iz*(h/2)*(Z*e3+S))/e3;
% plot(Z,tao3,'b.-','LineWidth',2,'MarkerSize',20)
% hold on;
% plot(Z,t3,'r.-','LineWidth',2,'MarkerSize',20)
% grid on;
% set(gca,'XMinorTick','on');
% xlabel('Axis Z (mm)','FontWeight','bold');
% ylabel('La contrainte (MPa)','FontWeight','bold');
% legend('Cas 2(sans semelle)','Cas 4(avec semelle)')
% grid minor;
%%%%%%%%%%%%%%%%%%%%%%SEMELLE
% 
% 
% 
% e1=1;
% e2=1;
% e3=1;
% b=500;
% h=400;
% w=160;
% x=60;
% S=9600;
% TY=10000;
% Iz1=b*e3^3+b*e3*h^2+e1*h^3/6+pi*e2*h^3/8;
% Iz2=4*(w*(x^3)/12+w*x*(h/2)^2);
% Iz=Iz1+Iz2;
% 
% syms Y;
% t1M = e1*h^2*(1/6*((2+3*pi*(e3*b)/(e2*h)+3*pi*S/(e2*h))/(2+pi*e1/e2))+1/8*(1-4*Y^2/h^2))*TY/Iz;
% vpa(t1M)
% Y=-200:20:200;
% t1M = e1*h^2*(1/6*((2+3*pi*(e3*b)/(e2*h)+3*pi*S/(e2*h))/(2+pi*e1/e2))+1/8*(1-4*Y.^2/h^2))*TY/Iz;
% plot(Y,t1M,'r-');
% 
% syms o;
% t2P = e2*h^2*(1/3*((3*e3/e1*b/h+3*S/e1/h-1)/(pi+2*e2/e1))+1/4*sin(o))*TY/Iz;
% vpa(t2P)
% o=0:1/18*pi:pi;
% t2P = e2*h^2*(1/3*((3*e3/e1*b/h+3*S/e1/h-1)/(pi+2*e2/e1))+1/4*sin(o))*TY/Iz;
% plot(o,t2P,'r-');
% 
% syms Z;
% t3=-TY/Iz*(h/2)*(Z*e3+S);
% vpa(t3)
% Z=-500:20:500;
% t3=-TY/Iz*(h/2)*(Z*e3+S);
% plot(Z,t3,'r-');
