%% Rounding errors in Singular Matrix 
%% QYQ 2017.7.3
tic;
clear;

%% set parameter

i=4;% the partial parameter Cij
j=100;% as the above
vt=0:0.01:10;
nf=200;
fmax=200;
fc=3;
alpha=-13/3;
%C=random('uniform',10^3,10^8,1000,1000);
C=powerlaw_cov(vt,nf,fmax,fc,alpha);
N=length(C);
x=ones(N,1);
r=rank(C);%trace of Matrix
C1=C;
ita=10^(-6);%10^(-9);%10^(-6);% relative error
dc=ita*C(i,j);% delta C
C1(i,j)=C(i,j)+dc;
dC=det(C);
iC=pinv(C);
k=cond(C);% condition number C in norm 2
%p=10;% p usually is less than 10 
dC1=det(C1);
iC1=pinv(C1);

%% numerical

lamda=-1/2*log(abs(det(C)))-1/2*x'*iC*x;
lamda2=-1/2*log(abs(det(C1)))-1/2*x'*iC1*x;
dlamda=lamda2-lamda;


%% Calculate the Matrix A

A=zeros(N,N);

for l=1:1:N
    for m=1:1:N
        A(l,m)=iC(l,i)*iC(j,m);
    end
end

ddlamda=abs((-1/2*iC(j,i)+1/2*x'*A*x)*dc);% analytic formula of delta lambda
rerr1=ddlamda/lamda;%relative error
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

%% Calculate Euclidian Norm for inv(C) 
% replace by norm(iC,'fro')
% a=0;
% for l=1:1:N
%     for m=1:1:N
%         a=iC(l,m)^2+a; 
%     end
%     b=sqrt(a); % the Euclid Norm
% end



%% error analysis
% if k >= 10^13
% [U,S,V]=svd(C);%SVD C
% s=pinv(S);% S is also singular
% merr=-2*eps*norm(V,'fro')*norm(U,'fro')*norm(s,'fro');% max error for inv(C)
% mmerr=(2*N+3)*eps*norm(U,'fro')^4*norm(x)^2*norm(s,'fro')^2;% error for multiple matrix times.
% error=merr+mmerr;% total error
% rerr2=error/C(i,j);
% 
% else 
%     merr=-p*k*eps*norm(iC,'fro');
%     mmerr=(p^2*k^2*eps^2+2*p*k*eps)*norm(iC,'fro')^2*norm(x)^2;
%     error=merr+mmerr;
%     rerr2=error/C(i,j);
% end
[U,S,V]=svd(C);
[U1,S1,V1]=svd(C1);
% Y=U'*x;
% s=S;
% s(find(s~=0))=1./S(find(S~=0));
dS=S1-S;
% delta=eye(N,N);
% 
% for l=1:1:N
%     delta(l,l)=dS(l,l)/S(l,l);
% end

error=2*N^2*ita;

toc;


