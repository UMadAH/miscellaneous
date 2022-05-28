%% PDE.Laplace Function with Circular boundary
clc
close all
%% global parameters, initial matrices
global h;
global N;
global a;
h = pi./N;
N = 200;
a = 20;
IC = 2;
Alpha = 10^(-4);
K_ij= zeros(2*N,2*N);
A_ij = zeros(2*N,2*N);
f_i = zeros(2*N,1);

%% define K_ij matrix£¨as t=s ,t/=s£©
for i= 1:2*N
    for j= 1:2*N
     if i==j
       K_ij(i,j) = -1 ./ pi .* log(a);
     else
       Hi=h*i;
       Hj=h*j;
       z = 0;
       z = sqrt(abs(a^2 * (cos(Hi)-cos(Hj))^2 + a^2 * (sin(Hi)-sin(Hj))^2));
       z = double(z);
       K_ij(i,j)= - 1 / (2*pi) * log(z /(4.0*(sin((Hi-Hj)/2))^2));
     end  
    end
end

%% define A_ij;R_ij matrices
for i=[1:2*N]
    for j=[1:2*N]
    r_ij=0;
rr_ij = 1 ./(2* N ).*(cos(N.*h*(i-j)));
    for m = 1:N-1
        r_ij = 1./m .* cos(m.*h*(i-j)) +r_ij;
        
    end
if i==j  
    A_ij(i,j) = 1/N*(r_ij+rr_ij)- 1./ N  .* log(a);
else
    A_ij(i,j) = 1/N*(r_ij+rr_ij)+ pi ./ N * K_ij(i,j);
end
    end
end 

%% define f_i matrix
for i = 1:2*N
    x1 = cos(h*(i-1));
    x2 = sin(h*(i-1));
    f_i(i,1)=(sin(a .*x1)+cos(a .*x2));
end

%% Solve linear equations (different iterations)
PSI = zeros(2*N,1);
I = eye(2*N);
aa =  1/norm(A_ij)^2 * rand(1,1);
if IC == 0
    PSI = inv(A_ij)*f_i;
elseif IC == 1
        for m = 1:100000
           PSI = PSI + (aa).* A_ij'* (f_i - A_ij * PSI);
        end
elseif IC == 2  
       PSI = inv(Alpha * I + A_ij' * A_ij)* A_ij' * f_i;
end            

%% solution&plot
AXEx = linspace(0, 2*pi, 2*N);
AXEy = linspace(0,a,2*N);
[AXE1, AXE2] = meshgrid(AXEx, AXEy); 
[AXE11, AXE22] = pol2cart(AXE1, AXE2); 
[AXE1, AXE2] = meshgrid(AXEx,[a]);
[AXEX, AXEY] = pol2cart(AXE1,AXE2);

for i =1:length(AXEy)
    for j =1:length(AXEx)
    CAN = (AXE11(i,j).*ones(1,length(AXEX))-AXEX);
    YYY = (AXE22(i,j).*ones(1,length(AXEY))-AXEY);
    U(i,j) = -sum(log(sqrt(sum(CAN.^2+YYY.^2,1))).*PSI'.*h/2,2)/pi;
    end
end
% surf(AXE11,AXE22,U);
% shading interp
% title('Default Comparison');
% axis equal
% hold on
% plot3(AXEX,AXEY,f_i);
figure
plot3(AXEX,AXEY,f_i)
hold on
plot3(AXE11(399,:),AXE22(399,:),U(399,:))
legend("f","U",'FontSize',12,'Location','Northwest')
title('Boundary Value Comparison');
axis equal

