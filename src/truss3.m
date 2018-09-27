% 3D Truss Analysis
% Input Data
NN=xlsread('truss3.xls','Sheet6','F3');
Coor=xlsread('truss3.xls','Sheet6','C9:F100');
for j=1:NN
    i=Coor(j,1);
    X(i)=Coor(j,2);
    Y(i)=Coor(j,3);
    Z(i)=Coor(j,4);
end

NE=xlsread('truss3.xls','Sheet6','F4');
Pro=xlsread('truss3.xls','Sheet6','H9:L100');
for j=1:NE
    i=Pro(j,1);
    E(i)=Pro(j,2);
    A(i)=Pro(j,3);
    Dir(i,1)=i;
    Dir(i,2)=Pro(j,4);
    Dir(i,3)=Pro(j,5);
end

NC=xlsread('truss3.xls','Sheet6','R3');
Cons=xlsread('truss3.xls','Sheet6','P9:Q100');

NF=xlsread('truss3.xls','Sheet6','Y3');
Forc=xlsread('truss3.xls','Sheet6','V9:Y100');


% Defination of Global Stiffness Matrix Element for each element
for i=1:NE
    L(i)=sqrt((X(Dir(i,3))- X(Dir(i,2)))^2+(Y(Dir(i,3))-Y(Dir(i,2)))^2+(Z(Dir(i,3))-Z(Dir(i,2)))^2);
    cosx(i)=(X(Dir(i,3))-X(Dir(i,2)))/L(i);
    cosy(i)=(Y(Dir(i,3))-Y(Dir(i,2)))/L(i);
    cosz(i)=(Z(Dir(i,3))-Z(Dir(i,2)))/L(i);
end

for i=1:NE
    S(i)=E(i)*A(i)/L(i);
    K11(:,:,i)=S(i)*[(cosx(i))^2 cosx(i)*cosy(i) cosx(i)*cosz(i);cosx(i)*cosy(i) (cosy(i))^2 cosy(i)*cosz(i);cosx(i)*cosz(i) cosy(i)*cosz(i) (cosz(i))^2];
    K12(:,:,i)=-K11(:,:,i);
    K21(:,:,i)=K12(:,:,i);
    K22(:,:,i)=K11(:,:,i);
end

% Defination of Structure Stiffness Matrix
K=zeros(3*NN,3*NN);
for n=1:NE
    i=Dir(n,2);
    j=Dir(n,3);
    K(3*i-2:3*i,3*i-2:3*i)=K11(:,:,n)+K(3*i-2:3*i,3*i-2:3*i);
    K(3*i-2:3*i,3*j-2:3*j)=K12(:,:,n);
    K(3*j-2:3*j,3*i-2:3*i)=K21(:,:,n);
    K(3*j-2:3*j,3*j-2:3*j)=K22(:,:,n)+K(3*j-2:3*j,3*j-2:3*j);
end

% Defination of Primary Nodal Forces
F=zeros(3*NN,1);
for i=1:NF
    f=3*Forc(i,1);
    F(f-2,1)=Forc(i,2);
    F(f-1,1)=Forc(i,3);
    F(f,1)=Forc(i,4);
end

% Elimination of rows and columns of K-matrix with respect to Supports
S=K;
for i=1:NC;
    r=3*Cons(i,1);
    if Cons(i,2)==0;
        S(r-2,:)=0; S(:,r-2)=0; S(r-2,r-2)=1;
        S(r-1,:)=0; S(:,r-1)=0; S(r-1,r-1)=1;
        S(r,:)=0; S(:,r)=0; S(r,r)=1;
    elseif Cons(i,2)==1;
        S(r-2,:)=0; S(:,r-2)=0; S(r-2,r-2)=1;
    elseif Cons(i,2)==2
        S(r-1,:)=0; S(:,r-1)=0; S(r-1,r-1)=1;
    elseif Cons(i,2)==3
        S(r,:)=0; S(:,r)=0; S(r,r)=1;
    end
end
% Solution of {F}={S}{D}
IS=pinv(S); 
I=IS;
for i=1:NC;
    r=3*Cons(i,1);
    if Cons(i,2)==0;
        I(r-2,:)=0; I(:,r-2)=0; I(r-2,r-2)=1;
        I(r-1,:)=0; I(:,r-1)=0; I(r-1,r-1)=1;
        I(r,:)=0; I(:,r)=0; I(r,r)=1;
    elseif Cons(i,2)==1;
        I(r-2,:)=0; I(:,r-2)=0; I(r-2,r-2)=1;
    elseif Cons(i,2)==2
        I(r-1,:)=0; I(:,r-1)=0; I(r-1,r-1)=1;
    elseif Cons(i,2)==3
        I(r,:)=0; I(:,r)=0; I(r,r)=1;
    end
end
d=I*F;

% Calculation of Nodal Forces Vector
W=K*d;

% Analysis Results
disp('3D-Truss Anlysis')
disp('Node Displacements')
disp(' ')
for i=1:NN
    fprintf('dx%g',i); fprintf('=%g\n',d(3*i-2))
    fprintf('dy%g',i); fprintf('=%g\n',d(3*i-1))
    fprintf('dz%g',i); fprintf('=%g\n',d(3*i))
end
disp('-------------------------------')


% Calculation of Elements Axial Force
for n=1:NE;
    i=Dir(n,2);
    j=Dir(n,3);
    PJ=K21(:,:,n)*[d(3*i-2);d(3*i-1);d(3*i)] + K22(:,:,n)*[d(3*j-2);d(3*j-1);d(3*j)];
    P(n)=PJ(1)*cosx(n)+PJ(2)*cosy(n)+PJ(3)*cosz(n);
    stress(n)=P(n)/A(n);
    strain(n)=stress(n)/E(n);
end
disp('Support Reactions:')
disp('')
for i=1:NC
    w=3*Cons(i,1);
    if Cons(i,2)==0
        fprintf('Rx%g\t',w/3); fprintf('=%g\n',W(w-2))
        fprintf('Ry%g\t',w/3); fprintf('=%g\n',W(w-1))
        fprintf('Rz%g\t',w/3); fprintf('=%g\n',W(w))
    elseif Cons(i,2)==1
        fprintf('Rx%g\t',w/3); fprintf('=%g\n',W(w-2))
    elseif Cons(i,2)==2
        fprintf('Ry%g\t',w/3); fprintf('=%g\n',W(w-1))
    else
        fprintf('Rz%g\t',w/3); fprintf('=%g\n',W(w-1))
 
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
 
 



