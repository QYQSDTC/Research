% Function to calculate the log likelihood ratio for particle swarm 
% optimazition (PSO) code developed by Prof. Soumya Mohanty.
% Adopt from LLR_PSOmpp.m, Average/Marginalize over pulsar phases. YW, Oct 24, 2016.

function [fitVal,varargout]=LLR_PSOav(xVec,inParams)
%function fitVal=LLR_PSO(xVec,inParams)

%rows: points
%columns: coordinates of a point
[nrows,ncols]=size(xVec);

%storage for fitness values
fitVal = zeros(nrows,1);

validPts = chkstdsrchrng(xVec);
% %Set fitness for invalid points to infty
fitVal(~validPts)=inf;
% if isempty(inParams)
%     %use default range of coordinates
%     %(only the numerical values below should be changed for different
%     %fitness functions)
%     xVec(validPts,:) = s2rscalar(xVec(validPts,:),-5,5);
% else
%     xVec(validPts,:) = s2rvector(xVec(validPts,:),inParams);
% end

%x=zeros(1,ncols);
% ===============================
%realCoord = zeros(size(xVec));
Np=inParams.Np;
realCoord = zeros(1,7+Np);
for lpc = 1:nrows
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);
        % 0<=x<=1 is convert to physical unit in 'avPhaseLLR'
            [ft,dummy]=avPhaseLLR(x,inParams);  % average over phases
            fitVal(lpc)= ft;
            realCoord(lpc,:)=dummy;
            %a=dummy1;
    end
end
% ===============================


%Return real coordinates if requested
if nargout > 1
    varargout{1}=realCoord;
end
% END of function