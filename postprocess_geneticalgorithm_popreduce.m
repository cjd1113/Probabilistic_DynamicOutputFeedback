


%This script performs post-processing of the GENETIC ALGORITHM probabilistic
%robust optimized controllers

%% THIS SCRIPT IS FOR THE PopReduce CASE %%

%Prepared by: Chris D'Angelo
%Date: September 5, 2018

clear all 
close all

warning('OFF','all')

%% LOAD CONTROLLER SIZE DATA USED FOR GA SYNTHESIS

cd '/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/Genetic Algorithms/FinalFunctionsforSuperComputing'

addpath('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/Genetic Algoritms/FinalFunctionsforSuperComputing')

initialpopulationloadingnopadding();


%%%%CHANGE THIS VARIABLE%%%%
%Declare the controller pair case that we are interested in.
Case = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Controller 1
controllersubstructure_index = 1; tolerance = 1e-3;

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


%Controller 2
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

%% LOAD PARTICLE GA RESULTS AND CONSTRUCT INTO DECENTRALIZED CONTROLLERS

%%%%%%ADJUST THIS PATH%%%%%%%%%
pathcase = sprintf('%d',Case);
pathname = fullfile('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/Genetic Algorithms/FinalFunctionsforSupercomputing/SupercomputerResults/PopReduce',pathcase);

cd(pathname)

filename = sprintf('GAOutputResults_Case%d2',Case);
load(filename,'-mat')

cd '/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/Genetic Algorithms/FinalFunctionsforSuperComputing'

Results

ControllerVariables = Results.Controllers;

ControllerVariables_1 = ControllerVariables(1:Controller1_length);
        
ControllerVariables_2 = ControllerVariables(Controller1_length+1:end);

[Ak1particle,Bk1particle,Ck1particle,Dk1particle] = ...
            controller_reconstruction(ControllerVariables_1,...
            numrealeigs1,nk1,mk1,rk1);
        
[Ak2particle,Bk2particle,Ck2particle,Dk2particle] = ...
            controller_reconstruction(ControllerVariables_2,...
            numrealeigs2,nk2,mk2,rk2);

Controller1Opt = struct('Ak',Ak1particle,'Bk',Bk1particle,'Ck',Ck1particle,'Dk',Dk1particle);
Controller2Opt = struct('Ak',Ak2particle,'Bk',Bk2particle,'Ck',Ck2particle,'Dk',Dk2particle);
        
%% FORMULATE THE CLOSED LOOP SYSTEM

%%CREATE BEAM #1

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

%%CREATE BEAM #2

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

%%FORM NOMINAL SYSTEM AND EVALUATE THE OPEN LOOP PERFORMANCE
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

coupledsystem = loopformulations_silent(beammodel_1,beammodel_2,Controller1Opt,Controller2Opt,coupling,'scenario_5');

Acoupled = coupledsystem.A;
Bcoupled = coupledsystem.B;
Ccoupled = coupledsystem.C;
Dcoupled = coupledsystem.D;


olsystem = loopformulations_silent(beammodel_1,beammodel_2,[],[],coupling,'scenario_0');
Aol = olsystem.A;
Bol = olsystem.B;
Col = olsystem.C;
Dol = olsystem.D;

coupledsystem_mu = loopformulations_silent(beammodel_1,beammodel_2,muController1,muController2,coupling,'scenario_5');
Aclmu = coupledsystem_mu.A;
Bclmu = coupledsystem_mu.B;
Cclmu = coupledsystem_mu.C;
Dclmu = coupledsystem_mu.D;


%% GENERATE MAXIMUM SINGULAR VALUE PLOTS 

w = logspace(-2,5,1000);
s = 1j.*w;
maxsv = zeros(1,length(s));
maxsvol = zeros(1,length(s));
maxsvmu = zeros(1,length(s));
for i = 1:length(s)
    
   Gprob(:,:,i) = Ccoupled/(s(i)*eye(size(Acoupled,1)) - Acoupled)*Bcoupled + Dcoupled;
   [~,S,~] = svd(Gprob(:,:,i));
   maxsv(i) = max(diag(S));
    
   Gol(:,:,i) = Col/(s(i)*eye(size(Aol,1)) - Aol)*Bol + Dol;
   [~,So,~] = svd(Gol(:,:,i));
   maxsvol(i) = max(diag(So));
   
   Gmu(:,:,i) = Cclmu/(s(i)*eye(size(Aclmu,1)) - Aclmu)*Bclmu + Dclmu;
   [~,Smu,~] = svd(Gmu(:,:,i));
   maxsvmu(i) = max(diag(Smu));
    
end

wHz = (1/(2*pi)).*w;

figure(1)
semilogx(wHz,20*log10(maxsvol),wHz,20*log10(maxsv),wHz,20*log10(maxsvmu),'linewidth',2)
leg = legend('Open Loop','Probabilistic Robust','$\mu$-Synthesized');
leg.FontSize = 14;
set(leg,'Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex','fontsize',14)
ylabel('Maximum Singular Value (dB)','interpreter','latex','fontsize',14)
title('Maximum Singular Value Plot - Open Loop, $\mu$-Synthesized, and Probabilistic Robust Controllers','interpreter','latex','fontsize',16)
set(gcf,'Units','Normalized','OuterPosition',[0 0.04 0.5 0.5])
    
cd '/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/Genetic Algorithms/FinalFunctionsforSuperComputing/SupercomputerResults/PopReduce/Figures'

figname = sprintf('GASolutionPopReduce_Case%d',Case);
saveas(gcf,figname,'png')

figure(2)
semilogx(wHz,20*log10(maxsv),wHz,20*log10(maxsvmu),'linewidth',2)
leg = legend('Probabilistic Robust','$\mu$-Synthesized');
leg.FontSize = 14;
set(leg,'Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex','fontsize',14)
ylabel('Maximum Singular Value (dB)','interpreter','latex','fontsize',14)
title('Comparison between $\mu$-Synthesized and Probabilistic Robust Controllers','interpreter','latex','fontsize',16)
set(gcf,'Units','Normalized','OuterPosition',[1 0.96 0.5 0.5])

figname2 = sprintf('GASolutionCompare_NoOL_Case%d',Case);    
saveas(gcf,figname2,'png')
