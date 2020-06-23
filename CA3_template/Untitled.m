%for n=1 case of Hagen Poiseulli
clear 
clc
n=5; %number of elements (only odd)
nn=n+2; %number of nodes
m=(nn+1)/2; %centre node
R=1; %radius 
h=R/n; %size of each element
r(nn)=0;r(1)=R;
b(m)=R/2;
for i=1:n 
    r(nn-i)=(i-1)*h;
end
for j=2:nn-1
    C(j)=-5e3*(r(j-1)^2-r(j)^2);
end
C(1)=0; C(nn)=0;
A=zeros(nn,nn);
A(1,1)=1/h; A(nn,nn)=-1/h;A(nn,nn-1)=1/h;
for i=2:nn-1
    r(1)=2*R;
    A(i,i-1)=(r(i-1))/h;
    A(i,i+1)=(r(i))/h;
    A(i,i)=(-A(i,i-1)-A(i,i+1));
end
X = zeros(nn,1);
Error_eval = ones(nn,1);
%% Start the Iterative method
iteration = 0;
q=100; %number of iterations
while max(Error_eval) > 0.001
    iteration = iteration + 1;
    Z = X;  % save current values to calculate error later
    for i = 1:nn
        j = 1:nn; % define an array of the coefficients' elements
        j(i) = [];  % eliminate the unknow's coefficient from the remaining coefficients
        Xtemp = X;  % copy the unknowns to a new variable
        Xtemp(i) = [];  % eliminate the unknown under question from the set of values
        X(i,1) = (C(i) - sum(A(i,j) * Xtemp)) / A(i,i);
    end
    Xsolution(:,iteration) = X;
    Error_eval = sqrt((X - Z).^2);
end
r(nn-1)=0;r(1)=R;
for k=1:(n-1)/2
    b(m-k)=b(m)+k*h;
    b(m+k)=b(m)-k*h;
end
b(1)=1;b(nn)=0;
%analytical result
u=2500*(1-b.^2);
plot(X,b)
hold on
plot(u,b)
legend('FVM value','Analytical value','location','NorthEastOutside')