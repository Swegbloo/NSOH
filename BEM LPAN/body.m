clc;
clear all;

load("excel.mat");

N = excel(:,1);
n = size(N);
n = n(:,1);
n2 = 1;

for i = 2:n
    if N(i)==1
        break;
    else
        n2 = N(i); %final value of n2 is the number of coordinates for any file of same format
    end
end

A = zeros(n2);
B = zeros(n2,1);

X = excel(1:n2,2);
Y = excel(1:n2,3);
Z = excel(1:n2,4);

t = 0;

CG = zeros(excel(n,1),3);

n3 = size(CG);
n3 = n3(1,1); %sus
v1 = zeros(n3,3);
v2 = zeros(n3,3);
n_vec = zeros(n3,3);
s = zeros(n3,1);

for i=1:n3
    CG(i,1) = (X(excel(n2+i,2))+X(excel(n2+i,3))+X(excel(n2+i,4))+X(excel(n2+i,5)))/4;
    CG(i,2) = (Y(excel(n2+i,2))+Y(excel(n2+i,3))+Y(excel(n2+i,4))+Y(excel(n2+i,5)))/4;
    CG(i,3) = (Z(excel(n2+i,2))+Z(excel(n2+i,3))+Z(excel(n2+i,4))+Z(excel(n2+i,5)))/4;
    v1(i,1) = X(excel(n2+i,2))-X(excel(n2+i,3));
    v2(i,1) = X(excel(n2+i,4))-X(excel(n2+i,3));
    v1(i,3) = Z(excel(n2+i,2))-Z(excel(n2+i,3));
    v2(i,3) = Z(excel(n2+i,4))-Z(excel(n2+i,3));
    v1(i,2) = Y(excel(n2+i,2))-Y(excel(n2+i,3));
    v2(i,2) = Y(excel(n2+i,4))-Y(excel(n2+i,3));
    abs_n = sqrt((v1(i,2)*v2(i,3)-v1(i,3)*v2(i,2))^2+(v1(i,3)*v2(i,1)-v1(i,1)*v2(i,3))^2+(v1(i,1)*v2(i,2)-v1(i,2)*v2(i,1)));
    n_vec(i,1) = (v1(i,2)*v2(i,3)-v1(i,3)*v2(i,2))/abs_n;
    n_vec(i,2) = (v1(i,3)*v2(i,1)-v1(i,1)*v2(i,3))/abs_n;
    n_vec(i,3) = (v1(i,1)*v2(i,2)-v1(i,2)*v2(i,1))/abs_n;
    s(i) = abs_n/2;
end

for i=1:n3
    
end

R_abs = zeros(n3);
R_vec = zeros(n3,n3,3);

for i = 1:n3
    for j = 1:n3
        R_abs(j,i) = sqrt((CG(j,1)-CG(i,1))^2+(CG(j,2)-CG(i,2))^2+(CG(j,3)-CG(i,3))^2);
        R_vec(j,i,1) = CG(j,1)-CG(i,1);
        R_vec(j,i,2) = CG(j,2)-CG(i,2);
        R_vec(j,i,3) = CG(j,3)-CG(i,3);
    end
end

A = get_inf_coef(R_vec,n_vec,s,R_abs);
B = get_B(R_vec,n_vec,s,R_abs);