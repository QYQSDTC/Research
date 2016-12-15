function [] = qwikperfeval(simDataDir,searchParamsFile)
% Performance evaluation of PSO based PTA analysis-quick and dirty version
% Same as PERFEVAL except that the number of PSO iterations is drastically
% reduced. Use this for debugging and testing code.


%-------------------------
%   QUICK PSO
%-------------------------

psoParams = struct('popSize', 5,...
                   'maxSteps', 10,...
                   'fminMaxFunEvals',10);
           
tic;

%%
%% Set the ranges of the search parameters
% The PSO search space is the same for all input files.
load(searchParamsFile,'xmaxmin');
disp(xmaxmin);
% Search Space dimensionality
nDim = length(xmaxmin(:,1));

%%
% Path to folder containing simulated data
% simDataDir = 'simData3';
% dir(simDataDir);
%%
% Path to folder for storing outputs
outDataDir = [simDataDir,filesep,'results'];
mkdir(outDataDir);
%%
% Basic ascii log file for storing run information
% logFile = [simDataDir,filesep,'RunLog.txt'];
% fidLog = fopen(logFile,'w');

% List of input data files (we will only process .mat files)
inputFiles = dir([simDataDir,filesep,'*.mat']);
%inputFiles = dir([simDataDir,filesep,'snr2*']);

% Number of input Data files 
nInputFiles = length(inputFiles);
%% 
% Number of independent PSO runs per file
nRuns = 4;
%%
% PSO configuration parameter structure
% psoParams = psoparamstruct(1,'default');
% psoParams.tuningVars.numPart = 40;
% psoParams.tuningVars.maxSteps = 2000;
% psoParams.tuningVars.inertiaDecayLaw = 'linear';
% psoParams.tuningVars.inertiaDecayParams = [0.9,0.4,2000,0.2];
% psoParams.topology.topoScheme = 'ring';
% psoParams.topology.topoParams = 3;
% psoParams.convergeScheme = 'maxSteps';
%------------------Best so far-----------------
% optStruct = optimset('fminsearch');
% optStruct.MaxFunEvals=200;
% psoParams = psoparamstruct(1,'default');
% psoParams.modScheme=struct('schemeName','nelderMead_gbest',...
%                                   'schemeParams',optStruct);
% psoParams.tuningVars.numPart = 40;
% psoParams.tuningVars.maxSteps = 2000;
% psoParams.convergeScheme = 'maxSteps';
% psoParams.tuningVars.inertiaDecayLaw = 'linear';
% psoParams.tuningVars.inertiaDecayParams = [0.9,0.4,psoParams.tuningVars.maxSteps,0.2];
% psoParams.topology.topoScheme = 'ring';
% psoParams.topology.topoParams = 3;
%------------------------------------
% PSO output parameter structure
% outP = struct('outFileNames',{{'log',''}},'graphics','off','status','off');

   
%% Process each input file
parfor lpFile = 1:nInputFiles
    
    inFileName = inputFiles(lpFile).name;
    %Output files carry the same name as input file
    outFileName = [outDataDir,filesep,inFileName];
    
    %Log info
    %fprintf(fidLog,'Starting on File %s.\n--------\n',inFileName);
    disp(['Starting on File ',inFileName]);
     
    %Load contents of input file
    inFileContents = load([simDataDir,filesep,inFileName]);
    %Data to be analyzed
    %data = inFileContents.timingResiduals;
    
    %Parameter structure for fitness function
    inParams = inFileContents.simParams;  
    % Add more info 
    yr = inFileContents.yr;
    timingResiduals = inFileContents.timingResiduals;
    inParams.s=timingResiduals;
    inParams.yr=yr;
    inParams.xmaxmin=xmaxmin;      
    Np=inParams.Np;
     
    % Target fitness function
    fHandle = @(x) LLR_PSOmpp(x,inParams);
    
    %Start independent PSO runs
    bestLocationVec = zeros(nRuns,nDim);
    bestLocRealC = zeros(nRuns,nDim+Np);
    bestFitValVec = zeros(nRuns,1);
    nFuncEvalsVec = zeros(nRuns,1);
    nIterVec = zeros(nRuns,1);
    wallClkTimeVec = zeros(nRuns,1);
    
    realC=zeros(1,nDim+Np);
    bestFitVal=0.0;
    
    %Log info
%     fprintf(fidLog,'\t fitness found/number of function evaluations \n ----\n');
%     disp('fitness found/number of function evaluations');
    for lprun = 1:nRuns
        
        tic;
        psoResults = ptapso(fHandle,nDim,psoParams);
        wallClkTime = toc;
        wallClkTimeVec(lprun) = wallClkTime/60;%min
        
        bestLocationVec(lprun,:)=psoResults.bestLocation;
        [~,bestLocRealC(lprun,:)] = LLR_PSOmpp(psoResults.bestLocation,inParams);
        bestFitValVec(lprun) = psoResults.bestSNR;
        nFuncEvalsVec(lprun) = psoResults.totalFuncEvals;
        nIterVec(lprun) = psoResults.totalSteps;
        
        %Log info
%         fprintf(fidLog,'\t %e/', bestFitValVec(lprun));
%         fprintf(fidLog,'%u', ...
%             nFuncEvalsVec(lprun));
        disp(['finishing run ',num2str(lprun)]);

    end
    
    [bestFitVal,bestFitIndx] = min(bestFitValVec);
%     fprintf(fidLog,'\n\t Best Run: %d \n',bestFitIndx);
    bestLocation = bestLocationVec(bestFitIndx,:);
    
    [~,realC]=LLR_PSOmpp(bestLocation,inParams);
    % disp('Best Coordinates (before fminsearch)');
    % disp(realC');
    % disp(['best fitness found (before fminsearch)', num2str(bestFitVal)]);
    % % Do fminsearch (Nelder-Mead)
    % bestLocation = fminsearch(fHandle,bestLocation);
    % [bestFitVal,realC]=LLR_PSO(bestLocation,inParams);
    % disp(['best fitness found (after fminsearch)', num2str(bestFitVal)]);
    % disp('Best Coordinates (after fminsearch)');
%     fprintf(fidLog,'\t Best location:\n');
%     fprintf(fidLog,'%f,',realC);
%     fprintf(fidLog,'\n--------\n');    
    
    %Make output file
    %--------------------
    outStruct = struct(...
                     'inputFile',[simDataDir,filesep,inFileName],...
                     'genHypothesis',inFileContents.genHypothesis,...
                     'id',inFileContents.id,...
                     'nRuns',nRuns,...
                     'bestLocVals',bestLocationVec,...
                     'bestLocRealCVals',bestLocRealC,...
                     'fitnessVals',bestFitValVec,...
                     'numFitEvals',nFuncEvalsVec,...
                     'numIter',nIterVec,...
                     'time2complete',wallClkTimeVec,...
                     'bestRun',bestFitIndx,...
                     'bestLocation',bestLocation,...
                     'bestRealLoc',realC...
                     );
     parsavelocal(outFileName,outStruct);
    
end
%End of run log info
% fclose(fidLog);

toc;
               
                   
