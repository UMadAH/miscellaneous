%2020/10/13;18:03;MWagon;BY MathisWang
clc
clear
close all
%X=xlsread('\X.xlsx','Sheet1','d6:k36');% 
X=xlsread('\X.xlsx','Sheet2','d3:k30')% 
Y=xlsread('\Y.xlsx','Sheet1','b2:i17');
[m,n]=size(Y);
[mou,n1]=size(X);

y_a=mean(Y');
y_a=y_a';

for i=1:m
    y_f=Y(i,1);
end
%for i=1:m
%    for j=1:n
%        y2(i,j)=y(i,j)/y_f(1);
%    end
%end   
for i=1:m
    for j=1:n
        y2(i,j)=Y(i,j)/y_a(i);
    end
end   

x_a=mean(X');
x_a=x_a';
for i=1:mou
    for j=1:n1
        x2(i,j)=X(i,j)/x_a(i);
    end
end   %X avg

for k=1:mou
%k=1,mou %choose specific target

    for i=1:m
        for j=1:n
            DELTA(i,j)=abs(x2(k,j)-y2(i,j));
        end
    end   %sequence(0,1)

    %minmax
    mmax=max(max(abs(DELTA)));
    mmin=min(min(abs(DELTA)));
    rho=0.5;
    %relative coef
    ksi=((mmin+rho*mmax)./(abs(DELTA)+rho*mmax));

    %degree of relativity
    ksi_column_num=size(ksi,2);
    r(:,k)=sum(ksi,2)/ksi_column_num;
end
    xlswrite('Z2.xlsx',r);
    %[rs,rind]=sort(r,'descend');