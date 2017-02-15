% calculate the best-fit timing residuals based on the estimated intrinsic
% and extrinsic parameters. Similar to 'bestfitResidualsfunc'
% May 27, 2014.  Yan Wang

function bfrsd=bfResidualsfunc(allParam,inParams)

%% transfer parameters from structure inParams
Np = inParams.Np;
N = inParams.N;
%s = inParams.s;
%sd = inParams.sd;
alphaP = inParams.alphaP;
deltaP = inParams.deltaP;
kp = inParams.kp;
yr = inParams.yr;
%xmaxmin = inParams.xmaxmin;

alpha=allParam(1);
delta=allParam(2);
omega=allParam(3);
phi0=allParam(4);
Amp=allParam(5);
iota=allParam(6);
thetaN=allParam(7);
%phiI=zeros(Np,1);
phiI=allParam(8:length(allParam));

Phi=omega*yr;  % Phi=omega*t, N by 1 matrix
bfrsd=zeros(Np,N);

% sky location of source in Cartesian coordinate
k=zeros(1,3);  % unit vector pointing from SSB to source
%for i=1:1:Ns
k(1)=cos(delta)*cos(alpha);
k(2)=cos(delta)*sin(alpha);
k(3)=sin(delta);
%end

for i=1:1:Np
   
    theta=acos( dot(k(1,:),kp(i,:)) );  % only one source here
    alphatilde=alpha-alphaP(i);
    
    Pp=-cos(deltaP(i))^2*(1-2*cos(alphatilde)^2+cos(alphatilde)^2*cos(delta)^2)...
        +sin(deltaP(i))^2*cos(delta)^2-0.5*sin(2*deltaP(i))*cos(alphatilde)*sin(2*delta);
    
    Pc=2*cos(deltaP(i))*sin(alphatilde)*(cos(deltaP(i))*cos(alphatilde)*sin(delta)...
        -sin(deltaP(i))*cos(delta));
    
    Fp=Pp/(1-cos(theta));
    Fc=Pc/(1-cos(theta));
    
    %disp(['in testMul','Pp',num2str(Pp),'Pc',num2str(Pc),'theta',num2str(theta)])
    
    A=2*Amp*sqrt( (1+cos(iota)^2)^2*(Fp*cos(2*thetaN)-Fc*sin(2*thetaN))^2 ...
        + 4*cos(iota)^2*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN))^2 ); % A is different for pulsars
    
    %tmp=-2*cos(iota)/(1+cos(iota)^2)*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN))/(Fp*cos(2*thetaN)-Fc*sin(2*thetaN));
    % solve psi, atan or atan2 ?
    %psi=atan(tmp);
    psi=atan2( -2*cos(iota)*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN)),(1+cos(iota)^2)*(Fp*cos(2*thetaN)-Fc*sin(2*thetaN)) );
    
    sin_ts=sin(phi0+psi+Phi);  % N by 1 matrix
    cos_ts=cos(phi0+psi+Phi);
    
    B=0.5*A*sin(phi0)*sin_ts;
    C=-0.5*A*cos(phi0)*sin_ts;
    D=0.5*A*sin(phi0)*cos_ts;
    E=-0.5*A*cos(phi0)*cos_ts;

    bfrsd(i,:)=(B-E)*cos(2*phiI(i))+(C+D)*sin(2*phiI(i))+(B+E);  % Eq 12
    %timingResiduals_test2(i,:)=A*sin(phi0-phiI(i))*sin(phi0+phiI(i)+psi+Phi);  % Eq 9

end

% disp('In bfResidualsfunc: ')
% for i=1:1:Np
%     cn=cos(2*phiI(i));
%     sn=sin(2*phiI(i));
%     disp([i,cn,sn]);
% end


% END of function