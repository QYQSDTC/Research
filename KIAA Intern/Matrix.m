%% Matrix Calculate
%% QYQ 2017.7.3

clear;

%% set parameter

%x=random('normal',1,10,1,3);
x=[1 1 1 1];
C=random('uniform',1,10,4,4);
C1=C+[0 0 0 0;0 0 0 0;0.1 0 0 0;0 0 0 0];%random('uniform',1,10,3,3);
dC=det(C);
iC=inv(C);
dC1=det(C1);
iC1=inv(C1);

%% numerical

lamda=-1/2*log(dC)-1/2*x*iC*x';
lamda2=-1/2*log(dC1)-1/2*x*iC1*x';
dlamda=lamda2-lamda;


%% analytic

%v=zeros(3,3);
A=zeros(4,4);
i=3;% the partial parameter Cij
j=1;% as the above
dc=C1(i,j)-C(i,j);% delta C  
%v(i,j)=1;

%% Calculate the Matrix A

for l=1:1:4
    for m=1:1:4
        A(l,m)=iC(l,i)*iC(j,m);
    end
end

ddlamda=(-1/2*iC(j,i)+1/2*x*A*x')*dc;
%err=-1/2*iC(j,i)*dc;
%err2=1/2*x*A*x'*dc;


% for s=1:1:3
%     l=s;% l represents the i in partial Cij
%     for r=1:1:3
%         m=r;% the j in partial Cij
%         dc=C1(l,m)-C(l,m);
%         for i=1:1:3
%             for j=1:1:3
%                 A(i,j)=iC(i,l)*iC(m,j);
%             end
%         end
%         ddlamda=(-1/2*iC(m,m)+1/2*x*A*x')*dc+ddlamda;
%     end
% end
