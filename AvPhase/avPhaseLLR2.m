% The log likelihood ratio function which is averaged/marginalized over pulsar 
% phases, this function should be implemented as efficient as possible, since 
% it would be called for a large amount of times; for single source.
% 10/19/2016, Yan Wang: adopt and modify from LogLikelihoodRatioMP5()
% 10/25/2016, YW, subtract the maximum value of b in each quadrant to avoid 
% infinity in the integration of the exponentials

function [LLR,varargout]=avPhaseLLR2(x,inParams)

%% transfer parameters from structure inParams
Np = inParams.Np;
%N = inParams.N;
s = inParams.s;
sd = inParams.sd;
alphaP = inParams.alphaP;
deltaP = inParams.deltaP;
kp = inParams.kp;
yr = inParams.yr;
xmaxmin = inParams.xmaxmin;

% transform x from standard coordinates [0,1] (requested by PSO) to physical 
% coordinates for 7 intrinsic parameters
alpha=x(1)*(xmaxmin(1,1)-xmaxmin(1,2))+xmaxmin(1,2);  % [0, 2*pi]
delta=x(2)*(xmaxmin(2,1)-xmaxmin(2,2))+xmaxmin(2,2);  % [-pi/2, pi/2]
omega=x(3)*(xmaxmin(3,1)-xmaxmin(3,2))+xmaxmin(3,2);
phi0=x(4)*(xmaxmin(4,1)-xmaxmin(4,2))+xmaxmin(4,2);
Amp=10^(x(5)*(xmaxmin(5,1)-xmaxmin(5,2))+xmaxmin(5,2));
iota=x(6)*(xmaxmin(6,1)-xmaxmin(6,2))+xmaxmin(6,2);
thetaN=x(7)*(xmaxmin(7,1)-xmaxmin(7,2))+xmaxmin(7,2);

%% calculate c for each pulsar
Phi=omega*yr;  % Phi=omega*t, N by 1 matrix
c=zeros(Np,4);  % Eq. 25
% sin_ts=zeros(N,1);
% cos_ts=zeros(N,1);
% B=zeros(N,1);
% C=zeros(N,1);
% D=zeros(N,1);
% E=zeros(N,1);

% sky location of source in Cartesian coordinate
k=zeros(1,3);  % unit vector pointing from SSB to source
%for i=1:1:Ns
k(1)=cos(delta)*cos(alpha);
k(2)=cos(delta)*sin(alpha);
k(3)=sin(delta);
%end

% coefficients used in the numerical integration for the marginalization 
b=zeros(Np,6);
%bn=zeros(1,6);  % normalized by max(abs(b))
%bs=zeros(1,6);  % b*sign
norm=zeros(4,1);  % normalization factor
%norm2=zeros(4,1);  % second normalization factor

phiI=zeros(Np,1);  % placeholder for the case of marginalization
LLR=0.0;  % log likelihood ratio
%LR=0.0;  % likelihood ratio
LRn=zeros(4,1);  % normalized LR in quadrant

%sign=zeros(4,6);
sign=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0; ...
      1.0, -1.0, 1.0, -1.0, 1.0, 1.0; ...
      1.0, -1.0, -1.0, 1.0, 1.0, 1.0; ...
      1.0, 1.0, -1.0, -1.0, 1.0, 1.0];

intup=[pi/2,pi,3*pi/2,2*pi];  % up limits of quadrants
intlow=[0,pi/2,pi,3*pi/2];  % low limits of quadrants

for i=1:1:Np
    
    theta=acos( dot(k,kp(i,:)) );
    [c(i,:),inn]= cfunc( alpha,delta,alphaP(i),deltaP(i),theta,...
                           Amp,omega,iota,thetaN,phi0,Phi,s(i,:),sd(i) );
    % c matrix will not be used in the marginalized likelihood
                       
    b(i,1) = inn(3)-0.5*inn(9);  % const
    b(i,2) = inn(1)-inn(6);  % cos(2p)
    b(i,3) = inn(2)-inn(8);  % sin(2p)
    b(i,4) = -inn(5);  % cos(2p)*sin(2p)
    b(i,5) = -0.5*inn(4);  % cos(2p)^2
    b(i,6) = -0.5*inn(7);  % sin(2p)^2

    % quadrant I, 
    for j=1:1:4
        bs=b(i,:).*sign(j,:);
        norm(j) = max(bs);
        bn = bs-norm(j);  % bn<=0
        fun=@(x)exp(sign(j,2)*bn(2).*cos(x)+sign(j,3)*bn(3).*sin(x) ...
            +sign(j,4)*bn(4).*sin(x).*cos(x)+sign(j,5)*bn(5).*cos(x).^2 ... 
            +sign(j,6)*bn(6).*sin(x).^2);  % x=2*pulsar_phase
        LRn(j) = integral(fun,intlow(j),intup(j)) * exp(sign(j,1)*bn(1));  % normalized LR
        
    end
    norm1=6*norm;
    N=max(norm1);
    norm2=norm1-N;  % norm2<=0
    
    M=0;
    for j=1:1:4
        M=M+exp(norm2(j))*LRn(j);
    end

    LLR = LLR + N + log(M);  % LLR is additive, LR is not.

%     disp('In LogLikelihoodRatioMP: ')
%     %disp([alpha,delta,omega,phi0,Amp,iota,thetaN,phiI']');
%     disp([i,cos(2*phiI(i)),sin(2*phiI(i))]);
    
end

LLR = -LLR;  % return -log likelihood ratio

if nargout > 1
    varargout{1}=[alpha,delta,omega,phi0,Amp,iota,thetaN,phiI'];  % in real coord.
end

% END of function