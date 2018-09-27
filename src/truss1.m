% Input Data
NN=xlsread('truss1.xls','Sheet1','F3');
Coor=xlsread('truss1.xls','Sheet1','C9:E100');
for j=1:NN
    i=Coor(j,1);
    X(i)=Coor(j,2);
    Y(i)=Coor(j,3);
end

NE=xlsread('truss1.xls','Sheet1','F4');
Pro=xlsread('truss1.xls','Sheet1','H9:L100');
for j=1:NE
    i=Pro(j,1);
    E(i)=Pro(j,2);
    A(i)=Pro(j,3);
    Dir(i,1)=i;
    Dir(i,2)=Pro(j,4);
    Dir(i,3)=Pro(j,5);
end

NC=xlsread('truss1.xls','Sheet1','R3');
Cons=xlsread('truss1.xls','Sheet1','P9:Q100');

NF=xlsread('truss1.xls','Sheet1','Y3');
Forc=xlsread('truss1.xls','Sheet1','V9:X100');

% Defination of Global Stiffness Matrix Element for each element
for i=1:NE
    L(i)=sqrt((X(Dir(i,3))- X(Dir(i,2)))^2+(Y(Dir(i,3))-Y(Dir(i,2)))^2);
    cos(i)=(X(Dir(i,3))-X(Dir(i,2)))/L(i);
    sin(i)=(Y(Dir(i,3))-Y(Dir(i,2)))/L(i);
end

for i=1:NE
    S(i)=E(i)*A(i)/L(i);
    K11(:,:,i)=S(i)*[(cos(i))^2 sin(i)*cos(i);sin(i)*cos(i) (sin(i))^2];
    K12(:,:,i)=-K11(:,:,i);
    K21(:,:,i)=K12(:,:,i);
    K22(:,:,i)=K11(:,:,i);
end

% Defination of Structure Stiffness Matrix
K=zeros(2*NN,2*NN);
for n=1:NE
    i=Dir(n,2);
    j=Dir(n,3);
    K(2*i-1:2*i,2*i-1:2*i)=K11(:,:,n)+K(2*i-1:2*i,2*i-1:2*i);
    K(2*i-1:2*i,2*j-1:2*j)=K12(:,:,n);
    K(2*j-1:2*j,2*i-1:2*i)=K21(:,:,n);
    K(2*j-1:2*j,2*j-1:2*j)=K22(:,:,n)+K(2*j-1:2*j,2*j-1:2*j);
end

% Defination of Primary Nodal Forces
F=zeros(2*NN,1);
for i=1:NF
    f=2*Forc(i,1);
    F(f-1,1)=Forc(i,2);
    F(f,1)=Forc(i,3);
end

% Elimination of rows and columns of K-matrix with respect to Supports
S=K;
for i=1:NC;
    r=2*Cons(i,1);
    if Cons(i,2)==0;
        S(r-1,:)=0; S(:,r-1)=0; S(r-1,r-1)=1;
        S(r,:)=0; S(:,r)=0; S(r,r)=1;
    elseif Cons(i,2)==1;
        S(r-1,:)=0; S(:,r-1)=0; S(r-1,r-1)=1;
    else
        S(r,:)=0; S(:,r)=0; S(r,r)=1;
    end
end

% Solution of {F}={S}{D}
I=inv(S);
d=I*F;

% Calculation of Nodal Forces Vector
W=K*d;

% Calculation of Elements Axial Force
for n=1:NE;
    i=Dir(n,2);
    j=Dir(n,3);
    PJ=K21(:,:,n)*[d(2*i-1);d(2*i)] + K22(:,:,n)*[d(2*j-1);d(2*j)];
    P(n)=PJ(1)*cos(n)+PJ(2)*sin(n);
    stress(n)=P(n)/A(n);
    strain(n)=stress(n)/E(n);
end

% Analysis Results
disp('2D-Truss Anlysis')
disp('Node Displacements')
disp(' ')
for i=1:NN
    fprintf('dx%g',i); fprintf('=%g\n',d(2*i-1))
    fprintf('dy%g',i); fprintf('=%g\n',d(2*i))
end
disp('-------------------------------')

disp('Support Reactions:')
disp('')
for i=1:NC
    w=2*Cons(i,1);
    if Cons(i,2)==0
        fprintf('Rx%g\t',w/2); fprintf('=%g\n',W(w-1))
        fprintf('Ry%g\t',w/2); fprintf('=%g\n',W(w))
    elseif Cons(i,2)==1
        fprintf('Rx%g\t',w/2); fprintf('=%g\n',W(w-1))
    else
        fprintf('Ry%g\t',w/2); fprintf('=%g\n',W(w))
    end
end

disp('--------------------------------')

disp('Elements Force:')
disp('')
for i=1:NE
    fprintf('P%g\t',i); fprintf('%g\n',P(i))
end
disp('--------------------------------')
disp('')
disp('Elements Stress:')   
for i=1:NE
    fprintf('Stress%g\t',i); fprintf('%g\n',stress(i))
end
disp('--------------------------------')
disp('')
disp('Elements Strain:')   
for i=1:NE
    fprintf('Strain%g\t',i); fprintf('%g\n',strain(i))
end
disp('--------------------------------')
disp('finished')    