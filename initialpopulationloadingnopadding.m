

%Controller loading, extraction, and initial population creation.

%For this version, controllers are NOT padded with zeros to bring them up to
%the largest controller dimension for either K1 or K2.  

%Case 1
%E \in [0.01E,2E] with D12 weighted by 1e-2

load('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/DK Iteration/Initial_Pop_GAInput/E_01_2/D12 1e-2/DKOutputScen1.mat')
load('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/DK Iteration/Initial_Pop_GAInput/E_01_2/D12 1e-2/DKOutputScen2.mat')



%Get all controllers from this case

Ak1_case1 = SolutionOutScen1.K.A;
Bk1_case1 = SolutionOutScen1.K.B;
Ck1_case1 = SolutionOutScen1.K.C;
Dk1_case1 = SolutionOutScen1.K.D;

K1_case1 = SolutionOutScen1.K;

Ak2_case1 = {};
Bk2_case1 = {};
Ck2_case1 = {};
Dk2_case1 = {};
K2_case1 = {};
for i = 1:length(SolutionOutScen2.infoc)
% for i = length(SolutionOutScen2.infoc):-1:1
    
    if SolutionOutScen2.infoc{i}.Bnd == Inf
        continue
    end
   Atemp = SolutionOutScen2.infoc{i}.K.A; 
   Btemp = SolutionOutScen2.infoc{i}.K.B;
   Ctemp = SolutionOutScen2.infoc{i}.K.C;
   Dtemp = SolutionOutScen2.infoc{i}.K.D;
   
   Ak2_case1{i} = Atemp;
   Bk2_case1{i} = Btemp;
   Ck2_case1{i} = Ctemp;
   Dk2_case1{i} = Dtemp;
   
   K2_case1{i} = SolutionOutScen2.infoc{i}.K;
   
end


%Case 2
%E \in [0.01E,2E] with D12 weighted by 1e-3

load('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/DK Iteration/Initial_Pop_GAInput/E_01_2/D12 1e-3/DKOutputScen1.mat')
load('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/DK Iteration/Initial_Pop_GAInput/E_01_2/D12 1e-3/DKOutputScen2.mat')

%Get all controllers from this case

Ak1_case2 = SolutionOutScen1.K.A;
Bk1_case2 = SolutionOutScen1.K.B;
Ck1_case2 = SolutionOutScen1.K.C;
Dk1_case2 = SolutionOutScen1.K.D;

K1_case2 = SolutionOutScen1.K;

Ak2_case2 = {};
Bk2_case2 = {};
Ck2_case2 = {};
Dk2_case2 = {};
K2_case2 = {};
for i = 1:length(SolutionOutScen2.infoc)
    
    if SolutionOutScen2.infoc{i}.Bnd == Inf
        continue
    end
   Atemp = SolutionOutScen2.infoc{i}.K.A; 
   Btemp = SolutionOutScen2.infoc{i}.K.B;
   Ctemp = SolutionOutScen2.infoc{i}.K.C;
   Dtemp = SolutionOutScen2.infoc{i}.K.D;
   
   Ak2_case2{i} = Atemp;
   Bk2_case2{i} = Btemp;
   Ck2_case2{i} = Ctemp;
   Dk2_case2{i} = Dtemp;
   
   K2_case2{i} = SolutionOutScen2.infoc{i}.K;
end


%Case 3
load('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/DK Iteration/Initial_Pop_GAInput/E_05_15/D12 1e-2/DKOutputScen1.mat')
load('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/DK Iteration/Initial_Pop_GAInput/E_05_15/D12 1e-2/DKOutputScen2.mat')

%Get all controllers from this case

Ak1_case3 = SolutionOutScen1.K.A;
Bk1_case3 = SolutionOutScen1.K.B;
Ck1_case3 = SolutionOutScen1.K.C;
Dk1_case3 = SolutionOutScen1.K.D;

K1_case3 = SolutionOutScen1.K;

Ak2_case3 = {};
Bk2_case3 = {};
Ck2_case3 = {};
Dk2_case3 = {};
K2_case3 = {};
for i = 1:length(SolutionOutScen2.infoc)
    
    if SolutionOutScen2.infoc{i}.Bnd == Inf
        continue
    end
   Atemp = SolutionOutScen2.infoc{i}.K.A; 
   Btemp = SolutionOutScen2.infoc{i}.K.B;
   Ctemp = SolutionOutScen2.infoc{i}.K.C;
   Dtemp = SolutionOutScen2.infoc{i}.K.D;
   
   Ak2_case3{i} = Atemp;
   Bk2_case3{i} = Btemp;
   Ck2_case3{i} = Ctemp;
   Dk2_case3{i} = Dtemp;
   
   K2_case3{i} = SolutionOutScen2.infoc{i}.K;
end


%Case 4
load('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/DK Iteration/Initial_Pop_GAInput/E_05_15/D12 1e-3/DKOutputScen1.mat')
load('/Users/christopherdangelo/Dropbox/Research/Probabilistic Decentralized Robust Active Structural Control/DK Iteration/Initial_Pop_GAInput/E_05_15/D12 1e-3/DKOutputScen2.mat')

%Get all controllers from this case

Ak1_case4 = SolutionOutScen1.K.A;
Bk1_case4 = SolutionOutScen1.K.B;
Ck1_case4 = SolutionOutScen1.K.C;
Dk1_case4 = SolutionOutScen1.K.D;

K1_case4 = SolutionOutScen1.K;

Ak2_case4 = {};
Bk2_case4 = {};
Ck2_case4 = {};
Dk2_case4 = {};
K2_case4 = {};
for i = 1:length(SolutionOutScen2.infoc)
    
    if SolutionOutScen2.infoc{i}.Bnd == Inf
        continue
    end
   Atemp = SolutionOutScen2.infoc{i}.K.A; 
   Btemp = SolutionOutScen2.infoc{i}.K.B;
   Ctemp = SolutionOutScen2.infoc{i}.K.C;
   Dtemp = SolutionOutScen2.infoc{i}.K.D;
   
   Ak2_case4{i} = Atemp;
   Bk2_case4{i} = Btemp;
   Ck2_case4{i} = Ctemp;
   Dk2_case4{i} = Dtemp;
   
   K2_case4{i} = SolutionOutScen2.infoc{i}.K;
   
end


ControllerPairing_case1 = cell(length(K2_case1),2);
for i = 1:length(K2_case1)
   
    ControllerPairing_case1{i,1} = K1_case1;
    ControllerPairing_case1{i,2} = K2_case1{i};
    
end


ControllerPairing_case2 = cell(length(K2_case2),2);
for i = 1:length(K2_case2)
   
    ControllerPairing_case2{i,1} = K1_case2;
    ControllerPairing_case2{i,2} = K2_case2{i};
    
end


ControllerPairing_case3 = cell(length(K2_case3),2);
for i = 1:length(K2_case3)
   
    ControllerPairing_case3{i,1} = K1_case3;
    ControllerPairing_case3{i,2} = K2_case3{i};    
    
end

ControllerPairing_case4 = cell(length(K2_case4),2);
for i = 1:length(K2_case4)
   
    ControllerPairing_case4{i,1} = K1_case4;
    ControllerPairing_case4{i,2} = K2_case4{i};
    
end

%%AGGREGATE THE INITIAL POPULATION
ControllerPairing_total = [ControllerPairing_case1;...
    ControllerPairing_case2;...
    ControllerPairing_case3;...
    ControllerPairing_case4];


return

%%VECTORIZE EACH CONTROLLER AND STORE IN A MATRIX
Kvectorized = [];
for i = 1:size(ControllerPairing_total,1)
    
   Ak1temp = ControllerPairing_total{i,1}.A;
   Bk1temp = ControllerPairing_total{i,1}.B;
   Ck1temp = ControllerPairing_total{i,1}.C;
   Dk1temp = ControllerPairing_total{i,1}.D;
    
   K1vec =  [Ak1temp(:);Bk1temp(:);Ck1temp(:);Dk1temp(:)];
   
   Ak2temp = ControllerPairing_total{i,2}.A;
   Bk2temp = ControllerPairing_total{i,2}.B;
   Ck2temp = ControllerPairing_total{i,2}.C;
   Dk2temp = ControllerPairing_total{i,2}.D;
    
   K2vec = [Ak2temp(:);Bk2temp(:);Ck2temp(:);Dk2temp(:)];
   
   Kstacked = [K1vec;K2vec];
   
   Kvectorized(:,i) = Kstacked;
    
end

InitialPopulation = Kvectorized;

