%MATLAB CODE
clear
fid=fopen('FTP75.txt','r');
[x,count]=fscanf(fid,'%f',[2,inf]);
[M,N]=size(x);
label=2; %label as speed
t1=x(1,1:N);
%plot(t1,x(label,1:N));
%peak value retrive, delete invalid amplitude
%delete non peak points
T(1)=x(2,1);n=1;t(1)=t1(1);
for i=2:(N-1)
if x(label,i)>x(label,i-1)&x(label,i)>=x(label,i+1)
n=n+1;T(n)=x(label,i);t(n)=t1(i);
elseif x(label,i)<x(label,i-1)&x(label,i)<=x(label,i+1)
n=n+1;T(n)=x(label,i);t(n)=t1(i);
else
end
end
n=n+1;T(n)=x(label,N);t(n)=t1(N);
plot(t1, x(label,1:N));
hold on
plot(t,T,'.','MarkerSize',12);
%extrapolate data, count value as points