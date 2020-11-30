
%% PARALLEL POST-PROCESSING AND STABILITY AND PERFORMANCE DEGRADATION FUNCTION CONSTRUCTION

%Prepared by: Chris D'Angelo
%Date: September 13, 2018

% Add path so that script knows where to find funcitons

addpath('/home/cjd66/ProbabilisticAnalysis/Functions')

%% Set up parallel pool on cluster
pc = parcluster('local')
parpool(pc,str2num(getenv('SLURM_CPUS_ON_NODE')))
poolobj = gcp;
%addAttachedFiles(poolobj,{' ',' '});

feature('numcores')
% 
spmd
    %Turn all warnings off on all workers to avoid memory leaks
    warning('off','all')
end

%Final result to analyze.  We have final results 1, 2, 3, 4, and 5 that we
%can analyze.
FinalResultIndex = 1;

%Set stretching factor vector.
sigmainctotal = [0 0.2 0.4 0.6 0.8 1 1.4 1.6 2];

%% CREATE BEAM #1

%%YOUNG'S MODULUS DEF'N %%
E = 200E9; %Young's modulus, Pa

%%REMAINING PARAMETER DEFINITIONS
rho = 7.800; %mass density
b = 5; %width, cm
h = .5; %height, cm
TotalLength = 100; %Beam length, cm
Ne = 10; %Number of finite elements

CtrlLoc = [10]; %Control input location, cm
DistLoc = 40; %Dist input location, cm
MeasLoc = [20]; %Meas input location, cm

Nc = round(Ne*(CtrlLoc./TotalLength),0);
Nd = round(Ne*(DistLoc./TotalLength),0);
Meas = round(Ne*(MeasLoc./TotalLength),0);

constraintloc = 1; %constrain beginning of beam, meaning the final two states will be used for coupling

%Generate the beam model
[Mcc,Kcc,KE,ME,Le,PD,PC,Mm] = bernoullibeamFEMfuncboundarydef(E,rho,b,h,Ne,TotalLength,Nc,Nd,Meas,constraintloc);

%Add damping
zeta_r = 0.02; %2 percent damping in each mode
[Ccc,V,D] = addDAMPING(Mcc,Kcc,zeta_r);

%Generate the generalized plant in physical coordinates
[A1, Bw_1, Bu_1, Cz_1, Cy_1, D111, D121, D211] = generalizedplant(Mcc,Ccc,Kcc,PD,PC,...
Mm,V,D,'nofilter','physical');

Bptemp1 = [zeros(size(A1,1)/2-2,2);eye(2)]; %
Bp_1 = [zeros(size(A1,1)/2,2);Mcc\Bptemp1];
% Cq_1 = [zeros(2,size(A1,1)-2), eye(2)];
Cq_1 = [zeros(2,size(A1,1)/2-2), eye(2), zeros(2,size(A1,2)/2)];

%% CREATE BEAM #2

%%YOUNG'S MODULUS DEF'N %%
E = 200E9; %Young's modulus, Pa

%%REMAINING PARAMETER DEFINITIONS
rho = 7.800; %mass density
b = 5; %width, cm
h = .5; %height, cm
TotalLength = 150; %Beam length, cm
Ne = 15; %Number of finite elements

CtrlLoc = 75; %Control input location, cm
DistLoc = 100; %Dist input location, cm
MeasLoc = 50; %Meas input location, cm

Nc = round(Ne*(CtrlLoc./TotalLength),0);
Nd = round(Ne*(DistLoc./TotalLength),0);
Meas = round(Ne*(MeasLoc./TotalLength),0);

constraintloc = 1; %constrain beginning of beam, meaning the final two states will be used for coupling

%Generate the beam model
[Mcc,Kcc,KE,ME,Le,PD,PC,Mm] = bernoullibeamFEMfuncboundarydef(E,rho,b,h,Ne,TotalLength,Nc,Nd,Meas,constraintloc);

%Add damping
zeta_r = 0.02; %2 percent damping in each mode
[Ccc,V,D] = addDAMPING(Mcc,Kcc,zeta_r);

%Generate the generalized plant in physical coordinates
[A2, Bw_2, Bu_2, Cz_2, Cy_2, D112, D122, D212] = generalizedplant(Mcc,Ccc,Kcc,PD,PC,...
Mm,V,D,'nofilter','physical');

%Now, we must define or identify, rather, the input and output influence
%matrices that will couple through the interconnection stiffness

Bptemp2 = [zeros(size(A2,1)/2-2,2);eye(2)]; %
Bp_2 = [zeros(size(A2,1)/2,2);Mcc\Bptemp2];
% Cq_2 = [zeros(2,size(A2,1)-2), eye(2)];
Cq_2 = [zeros(2,size(A2,1)/2-2), eye(2), zeros(2,size(A2,2)/2)];

%% FORM NOMINAL SYSTEM AND EVALUATE THE OPEN LOOP RESPONSE 

%Form structures and create nominal stiffness element, create the open loop
%system, and calculate the open-loop infinity norm.
beammodel_1 = struct('A',A1,'B1',Bw_1,'B2',Bu_1,'B3',Bp_1,...
    'C1',Cz_1,'C2',Cy_1,'C3',Cq_1,'D11',D111,'D12',D121,'D21',D211);
beammodel_2 = struct('A',A2,'B1',Bw_2,'B2',Bu_2,'B3',Bp_2,...
    'C1',Cz_2,'C2',Cy_2,'C3',Cq_2,'D11',D112,'D12',D122,'D21',D212);

%Generate a nominal stiffness
[K11nom,K12nom,K21nom,K22nom] = NomStiffness(b,h,10,E);
%Assign stiffness parameter global variables
K_width = b;
K_height = h;
K_length = 10;
K_Enom = E;

%Create interconnection stiffness structure for passing these variables on
%to optimization
Interconnection = struct('K_width',K_width,'K_height',K_height,'K_length',K_length,'K_Enom',K_Enom);

K11 = -K11nom; K12 = -K12nom; K21 = -K21nom; K22 = -K22nom;
coupling = struct('K11',K11,'K12',K12,'K21',K21,'K22',K22);

coupledsystem = loopformulations_silent(beammodel_1,beammodel_2,[],[],coupling,'scenario_0');

A_olnom = coupledsystem.A;
B_olnom = coupledsystem.B;
C_olnom = coupledsystem.C;
D_olnom = coupledsystem.D;

OLsystem = ss(A_olnom,B_olnom,C_olnom,D_olnom);

OLInfNorm = norm(OLsystem,inf);
OLInfNormdB = 20*log10(OLInfNorm);

InfNormReduction = -6;

CLInfNorm_desired = OLInfNormdB + InfNormReduction;

%% DESIRED PERFORMANCE LEVEL AND GA FITNESS FUNCTION PARAMETER SETTING

gamma_desired = 10^(CLInfNorm_desired/20); %Desired closed loop infinity norm performance

%Genetic Algorithm Cost Function Global Parameter Assignments
alpha = 60; %weighting for probability of stability
Beta = 40; %weighting for probability of performance
%% THIS WILL BE CHANGED FOR RIGOROUS ANALYSIS
ep = 0.01;
d = 0.05;
Nmce = ceil(1/(2*ep^2)*log(2/d));
% Nmce = 2;
% Nmce = 1153; %Number of monte carlo evaluations
%% THIS WILL BE CHANGED FOR DEGRADATION FUNCTION CONSTRUCTION
sigma = 0.4*E; %variance on the elastic modulus in the interconnection
% sigma = sigma*sigmainc;
%Structure with performance and uncertainty-related variables, which we will pass into optimization
PerformanceVariables = struct('gamma_desired',gamma_desired,'alpha',alpha,'Beta',Beta,'Nmce',Nmce,'sigma',sigma);

%% DECLARATION OF SYSTEM VARIABLE SIZES FOR FITNESS FUNCTION EFFICIENCY

%We need to get the sizes of certain system matrices, and declare these as
%global variables so that the size(*,*) function is not called during every
%monte carlo evaluation during the genetic algorithm.  This introduces
%unneccesary overhead.

%%Dynamic system
n_1 = size(A1,1); n_2 = size(A2,1);
m1_1 = size(Bw_1,2); m1_2 = size(Bw_2,2);
m2_1 = size(Bu_1,2); m2_2 = size(Bu_2,2);
r1_1 = size(Cz_1,1); r1_2 = size(Cz_2,1);
r2_1 = size(Cy_1,1); r2_2 = size(Cy_2,1);
%Zero out feedforward terms
D121 = 0*D121;
D211 = 0*D211;
D122 = 0*D122;
D212 = 0*D212;

%Create system size structures
system1sizes = struct('n_1',n_1,'m1_1',m1_1,'m2_1',m2_1,'r1_1',r1_1,'r2_1',r2_1);
system2sizes = struct('n_2',n_2,'m1_2',m1_2,'m2_2',m2_2,'r1_2',r1_2,'r2_2',r2_2);

%Assign zeroed-out feedforward matrices to beammodel structures.
beammodel_1.D12 = D121;
beammodel_1.D21 = D211;
beammodel_2.D12 = D122;
beammodel_2.D21 = D212;

WorkingDirectory = pwd;

%% LOADING CONTROLLER PAIRS FROM PROBABILISTIC ROBUST SYNTHESIS USING GENETIC ALGS OR PARTICLE SWARM OPTIMIZATION

%Note that we need this in the same directory, or have to cd to it.
% probabilisiticrobustcontrollerloading();
probabilisticrobustcontrollerloadingsupercomputer();

%Results are saved in a cell array that looks like FinalResults{1:5,1}
%The first four have the following structures:
%     Controllers: [1×1776 double]
%            FVAL: -99.8265
%        EXITFLAG: 1
%          OUTPUT: [1×1 struct]
%      POPULATION: [45×1776 double]
%          SCORES: [45×1 double]
%as these were found using genetic algorithm-based optimization.  The final
%has the structure:
%     Controllers: [1×1392 double]
%            FVAL: -74.1023
%        EXITFLAG: 1
%as it was found using particle swarm optimization

%Final results coincide with the following cases:
%FinalResults{1,1} - Case 7
%FinalResults{2,1} - Case 2
%FinalResults{3,1} - Case 4
%FinalResults{4,1} - Case 7
%FinalResults{5,1} - Case 2

%Note: we are going to have to use information derived from the initial
%seed controller for reconstructing our probabilistic robust solution.

ProbRobustControllerVector = FinalResults{FinalResultIndex,1}.Controllers;

%Now cd back into previous working directory
cd(WorkingDirectory)

%% LOADING CONTROLLER PAIRS FROM MU SYNTHESIS
initialpopulationloadingsupercomputer();
% initialpopulationloadingnopadding();

%Controller indicies and imaginary tolerance
%%%%CHANGE THIS VARIABLE%%%%
%Declare the controller pair case that we are interested in.
if FinalResultIndex == 1 || FinalResultIndex == 4 
    Case = 7;
elseif FinalResultIndex == 2 || FinalResultIndex == 5
    Case = 2;
elseif FinalResultIndex == 3
    Case = 4;
else
    disp('Invalid FinalResultIndex selection')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

controllersubstructure_index = 1; tolerance = 1e-3;

%%CONTROLLER #1

[Avectors_1,Bvectors_1,Cvectors_1,Dvectors_1,nk1,mk1,rk1,numrealeigs1,realeigindex] =...
    controllerdiagonalvectorization(ControllerPairing_total{Case,...
    controllersubstructure_index},tolerance);

[ControllerParameters_1,zeroimaginarypartindices] = ...
    controllervectorcompression(Avectors_1,Bvectors_1,Cvectors_1,Dvectors_1,numrealeigs1);

ControllerParameters_1(zeroimaginarypartindices) = [];

Controller1_length = length(ControllerParameters_1);

%%Extract mu synthesized controllers for comparison
Ak1mu = ControllerPairing_total{Case,controllersubstructure_index}.A;
Bk1mu = ControllerPairing_total{Case,controllersubstructure_index}.B;
Ck1mu = ControllerPairing_total{Case,controllersubstructure_index}.C;
Dk1mu = ControllerPairing_total{Case,controllersubstructure_index}.D;

muController1 = struct('Ak',Ak1mu,'Bk',Bk1mu,'Ck',Ck1mu,'Dk',Dk1mu);


%%CONTROLLER 2
controllersubstructure_index = 2; tolerance = 1e-3;

[Avectors_2,Bvectors_2,Cvectors_2,Dvectors_2,nk2,mk2,rk2,numrealeigs2,realeigindex] =...
    controllerdiagonalvectorization(ControllerPairing_total{Case,...
    controllersubstructure_index},tolerance);

[ControllerParameters_2,zeroimaginarypartindices] = ...
    controllervectorcompression(Avectors_2,Bvectors_2,Cvectors_2,Dvectors_2,numrealeigs2);

ControllerParameters_2(zeroimaginarypartindices) = [];

Controller2_length = length(ControllerParameters_2);

ControllerVariables = [ControllerParameters_1;ControllerParameters_2];

ControllerVariables_1 = ControllerVariables(1:Controller1_length);
        
ControllerVariables_2 = ControllerVariables(Controller1_length+1:end);

[Ak1,Bk1,Ck1,Dk1] = ...
            controller_reconstruction(ControllerVariables_1,...
            numrealeigs1,nk1,mk1,rk1);
        
[Ak2,Bk2,Ck2,Dk2] = ...
            controller_reconstruction(ControllerVariables_2,...
            numrealeigs2,nk2,mk2,rk2);

%%Extract mu synthesized controllers for comparison
Ak2mu = ControllerPairing_total{Case,controllersubstructure_index}.A;
Bk2mu = ControllerPairing_total{Case,controllersubstructure_index}.B;
Ck2mu = ControllerPairing_total{Case,controllersubstructure_index}.C;
Dk2mu = ControllerPairing_total{Case,controllersubstructure_index}.D;

muController2 = struct('Ak',Ak2mu,'Bk',Bk2mu,'Ck',Ck2mu,'Dk',Dk2mu);


%% Probabilistic robust controller reconstruction
%Reconstruct the probabilistically-robust controllers from the loaded
%probabilistically robust parameter vector and information about the seeds
%that we used during optimization.

ControllerVariables_1 = ProbRobustControllerVector(1:Controller1_length);
ControllerVariables_2 = ProbRobustControllerVector(Controller1_length+1:end);

[Ak1prob,Bk1prob,Ck1prob,Dk1prob] = ...
            controller_reconstruction(ControllerVariables_1,...
            numrealeigs1,nk1,mk1,rk1);
        
[Ak2prob,Bk2prob,Ck2prob,Dk2prob] = ...
            controller_reconstruction(ControllerVariables_2,...
            numrealeigs2,nk2,mk2,rk2);

Controller1Opt = struct('Ak',Ak1prob,'Bk',Bk1prob,'Ck',Ck1prob,'Dk',Dk1prob);
Controller2Opt = struct('Ak',Ak2prob,'Bk',Bk2prob,'Ck',Ck2prob,'Dk',Dk2prob);
        

E_nom = K_Enom;
L = K_length;
I = 1/12*K_width*K_height^3;

Adim = n_1+n_2+nk1+nk2;


%% FIRST WE WILL ANALYZE THE ROBUSTNESS AND GATHER DATA FOR THE MU-SYNTHESIZED SOLUTIONS

%Activate the mu synthesized solutions
Ak1 = Ak1mu; Bk1 = Bk1mu; Ck1 = Ck1mu; Dk1 = Dk1mu;
Ak2 = Ak2mu; Bk2 = Bk2mu; Ck2 = Ck2mu; Dk2 = Dk2mu;

for ii = 1:length(sigmainctotal)
    
    %Initialize counters
    CountStable = 0;
    CountPerf = 0;
    
    sigma = 0.4*Enom; %variance on the elastic modulus in the interconnection
    sigma = sigma*sigmainctotal(ii);
        
    %Memory preallocation
    % HinfNormData = zeros(Nmce,2);
    HinfNorm = zeros(Nmce,1);
    HinfNormfreq = zeros(Nmce,1);
    tic
    parfor i = 1:Nmce
        
        %Random interconnection stiffness element generation
        
        E = E_nom + sigma*randn;
        Ke = (E*I)/(L^3) * [12 6*L -12 6*L;...
            6*L 4*L^2 -6*L 2*L^2;...
            -12 -6*L 12 -6*L;...
            6*L 2*L^2 -6*L 4*L^2];
        
        K11 = -Ke(1:2,1:2);
        K12 = -Ke(1:2,3:4);
        K21 = -Ke(3:4,1:2);
        K22 = -Ke(3:4,3:4);
        
        %Formulate the closed loop system, which is "scenario 5" in the
        %loopformulations and loopformulations_silent functions.  Feedforward
        %terms are zeroed.
        
        
        A = [A1+Bu_1*Dk1*Cy_1+Bp_1*K11*Cq_1, Bp_1*K12*Cq_2, Bu_1*Ck1, zeros(n_1,size(Ak2,2));...
            Bp_2*K21*Cq_1, A2+Bu_2*Dk2*Cy_2+Bp_2*K22*Cq_2, zeros(n_2,size(Ck1,2)), Bu_2*Ck2;...
            Bk1*Cy_1, zeros(size(Bk1,1),n_2), Ak1, zeros(size(Ak1,1),size(Ak2,2));...
            zeros(size(Ak2,1),n_1), Bk2*Cy_2, zeros(size(Bk2,1),size(Ak1,2)), Ak2];
        
        B = [Bw_1+Bu_1*Dk1*D211, zeros(n_1,m1_2);...
            zeros(n_2,m1_1), Bw_2+Bu_2*Dk2*D212;...
            Bk1*D211, zeros(size(Bk1,1),m1_2);...
            zeros(size(Bk2,1),m1_1), Bk2*D212];
        
        C = [Cz_1+D121*Dk1*Cy_1, zeros(r1_1,n_2), D121*Ck1, zeros(r1_1,size(Ck2,2));...
            zeros(r1_2,n_1), Cz_2+D122*Dk2*Cy_2, zeros(r1_2,size(Ck1,2)), D122*Ck2];
        
        D = [D121*Dk1*D211, zeros(r1_1,m1_2);...
            zeros(r1_2,m1_1), D122*Dk2*D212];
        
        %Logical test for stability
        StabilityTest = all(real(eig(A))<0);
        CountStable = CountStable+StabilityTest;
        
        if StabilityTest == 1
            %         tic
            system = ss(A,B,C,D);
            [hinfnormval,infnormfreq] = norm(system,inf);
            HinfNorm(i) = hinfnormval;
            HinfNormfreq(i) = infnormfreq;
            %         toc
        else
            HinfNorm(i) = inf;
            HinfNormfreq(i) = inf;
        end
        
        %Form system and evaluate closed loop infinity norm
        
        %Sort natural frequencies into ascending real / imaginary parts
        
        if StabilityTest == 1
            
            clsystem_infnormviolated = 'false';
            cleigs = eig(A);
            
            %Only evaluate infinity norm at lightly damped modes
            zeta = abs(real(cleigs)./abs(cleigs));
            zetatol = 1e-4;
            zeta_crit = find(zeta>=1-zetatol);
            
            clevalfreqs = abs(cleigs);
            clevalfreqs(zeta_crit) = [];
            clevalfreqs_sorted = sort(clevalfreqs);
            clevalfreqs_sorttruncate = clevalfreqs_sorted(1:2:end/2); %Only evaluate up to middle of structure's bandwidth
            
            omega = clevalfreqs_sorttruncate;
            s = 1j.*omega;
            for k = 1:length(omega)
                G = C/(s(k)*eye(Adim) - A)*B+D;
                SingCheck = G'*G-(gamma_desired^2)*eye(size(B,2));
                if any(eig(SingCheck)>=0)
                    clsystem_infnormviolated = 'true';
                    break
                end
            end
            
            
            if strcmp(clsystem_infnormviolated,'false')==1
                CountPerf = CountPerf+1;
            end
            
        end
        
        
    end
    Elapsed_musyn = toc
    
    sigmainctotal(ii)
    
    Pstability = 1/Nmce*CountStable
    PPerf = 1/Nmce*CountPerf
    Cost = alpha*Pstability + Beta*PPerf
    
    HinfNormData = [HinfNorm HinfNormfreq];
    
    musyndata.PStability = Pstability;
    musyndata.PPerf = PPerf;
    musyndata.Cost = Cost;
    musyndata.sigma = sigma;
    musyndata.Nmce = Nmce;
    musyndata.HinfNormData = HinfNormData;
    musyndata.sigmainc = sigmainctotal(ii);
    
    %Save results
    
    dataname = sprintf('musyncontrollerresults_ProbRobCase%d_%3gEsigma.mat',FinalResultIndex,sigmainctotal(ii));
    save(dataname,'musyndata')

end

%clear results from mu syn data
clear Pstability PPerf Cost

%% Analysis of the probabilistic robust solution
%Activate the probabilistic robust solutions
Ak1 = Controller1Opt.Ak; Bk1 = Controller1Opt.Bk; Ck1 = Controller1Opt.Ck; Dk1 = Controller1Opt.Dk;
Ak2 = Controller2Opt.Ak; Bk2 = Controller2Opt.Bk; Ck2 = Controller2Opt.Ck; Dk2 = Controller2Opt.Dk;

for ii = 1:length(sigmainctotal)
    
    %Initialize counters
    CountStable = 0;
    CountPerf = 0;
    
    sigma = 0.4*Enom; %variance on the elastic modulus in the interconnection
    sigma = sigma*sigmainctotal(ii);
    
    %Memory preallocation
    % HinfNormData = zeros(Nmce,2);
    HinfNorm = zeros(Nmce,1);
    HinfNormfreq = zeros(Nmce,1);
    tic
    parfor i = 1:Nmce
        
        %Random interconnection stiffness element generation
        
        E = E_nom + sigma*randn;
        Ke = (E*I)/(L^3) * [12 6*L -12 6*L;...
            6*L 4*L^2 -6*L 2*L^2;...
            -12 -6*L 12 -6*L;...
            6*L 2*L^2 -6*L 4*L^2];
        
        K11 = -Ke(1:2,1:2);
        K12 = -Ke(1:2,3:4);
        K21 = -Ke(3:4,1:2);
        K22 = -Ke(3:4,3:4);
        
        %Formulate the closed loop system, which is "scenario 5" in the
        %loopformulations and loopformulations_silent functions.  Feedforward
        %terms are zeroed.
        
        
        A = [A1+Bu_1*Dk1*Cy_1+Bp_1*K11*Cq_1, Bp_1*K12*Cq_2, Bu_1*Ck1, zeros(n_1,size(Ak2,2));...
            Bp_2*K21*Cq_1, A2+Bu_2*Dk2*Cy_2+Bp_2*K22*Cq_2, zeros(n_2,size(Ck1,2)), Bu_2*Ck2;...
            Bk1*Cy_1, zeros(size(Bk1,1),n_2), Ak1, zeros(size(Ak1,1),size(Ak2,2));...
            zeros(size(Ak2,1),n_1), Bk2*Cy_2, zeros(size(Bk2,1),size(Ak1,2)), Ak2];
        
        B = [Bw_1+Bu_1*Dk1*D211, zeros(n_1,m1_2);...
            zeros(n_2,m1_1), Bw_2+Bu_2*Dk2*D212;...
            Bk1*D211, zeros(size(Bk1,1),m1_2);...
            zeros(size(Bk2,1),m1_1), Bk2*D212];
        
        C = [Cz_1+D121*Dk1*Cy_1, zeros(r1_1,n_2), D121*Ck1, zeros(r1_1,size(Ck2,2));...
            zeros(r1_2,n_1), Cz_2+D122*Dk2*Cy_2, zeros(r1_2,size(Ck1,2)), D122*Ck2];
        
        D = [D121*Dk1*D211, zeros(r1_1,m1_2);...
            zeros(r1_2,m1_1), D122*Dk2*D212];
        
        %Logical test for stability
        StabilityTest = all(real(eig(A))<0);
        CountStable = CountStable+StabilityTest;
        
        %Form system and evaluate closed loop infinity norm
        if StabilityTest == 1
            %         tic
            system = ss(A,B,C,D);
            [hinfnormval,infnormfreq] = norm(system,inf);
            HinfNorm(i) = hinfnormval;
            HinfNormfreq(i) = infnormfreq;
            %         toc
        else
            HinfNorm(i) = inf;
            HinfNormfreq(i) = inf;
        end
        %Sort natural frequencies into ascending real / imaginary parts
        
        if StabilityTest == 1
            
            clsystem_infnormviolated = 'false';
            cleigs = eig(A);
            
            %Only evaluate infinity norm at lightly damped modes
            zeta = abs(real(cleigs)./abs(cleigs));
            zetatol = 1e-4;
            zeta_crit = find(zeta>=1-zetatol);
            
            clevalfreqs = abs(cleigs);
            clevalfreqs(zeta_crit) = [];
            clevalfreqs_sorted = sort(clevalfreqs);
            clevalfreqs_sorttruncate = clevalfreqs_sorted(1:2:end/2); %Only evaluate up to middle of structure's bandwidth
            
            omega = clevalfreqs_sorttruncate;
            s = 1j.*omega;
            for k = 1:length(omega)
                G = C/(s(k)*eye(Adim) - A)*B+D;
                SingCheck = G'*G-(gamma_desired^2)*eye(size(B,2));
                if any(eig(SingCheck)>=0)
                    clsystem_infnormviolated = 'true';
                    break
                end
            end
            
            
            if strcmp(clsystem_infnormviolated,'false')==1
                CountPerf = CountPerf+1;
            end
            
        end
        
    end
    Elapsed_probrob = toc
    
    sigmainctotal(ii)
    
    Pstability = 1/Nmce*CountStable
    PPerf = 1/Nmce*CountPerf
    Cost = alpha*Pstability + Beta*PPerf
    
    HinfNormData = [HinfNorm HinfNormfreq];
    
    %Save results
    probrobdata.PStability = Pstability;
    probrobdata.PPerf = PPerf;
    probrobdata.Cost = Cost;
    probrobdata.sigma = sigma;
    probrobdata.Nmce = Nmce;
    probrobdata.HinfNormData = HinfNormData;
    probrobdata.sigmainc = sigmainctotal(ii);
    
    %Save results
    
    dataname = sprintf('probabilisticrobust_ProbRobCase%d_%3gEsigma.mat',FinalResultIndex,sigmainctotal(ii));
    save(dataname,'probrobdata')

        
end
