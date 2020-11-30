

%%LOOP AT A TIME TRANSFORMATION CODE

%Given two systems, modeled as generalized regulators with interconnection
%terms explicitly identified in input/output form, this code will formulate
%any/all system cases.  The scenario must be defined for this function.

%%THIS IS THE SILENT VERSION
%Details: The silent version does not display anything to the command
%prompt.  This is especially useful for when we are performing Monte Carlo
%evaluations and must formulate this system multiple times.

%Prepared by: Chris D'Angelo
%Date: May 31, 2018


function system = loopformulations_silent(system1,system2,controller1,controller2,coupling,scenario)


%Unpack the system structures

%% SYSTEM 1

A1 = system1.A;
Bw1 = system1.B1;
Bu1 = system1.B2;
Bp1 = system1.B3;
Cz1 = system1.C1;
Cy1 = system1.C2;
Cq1 = system1.C3;
D111 = system1.D11;
D121 = system1.D12;
D211 = system1.D21;


%% SYSTEM 2

A2 = system2.A;
Bw2 = system2.B1;
Bu2 = system2.B2;
Bp2 = system2.B3;
Cz2 = system2.C1;
Cy2 = system2.C2;
Cq2 = system2.C3;
D112 = system2.D11;
D122 = system2.D12;
D212 = system2.D21;


%% COUPLING COEFFICIENTS

K11 = coupling.K11;
K12 = coupling.K12;
K21 = coupling.K21;
K22 = coupling.K22;

%% CONTROLLERS

%Controller 1

if isempty(controller1)
    
    Ak1=0; Bk1=0; Ck1=0; Dk1=0;
    
else

    Ak1 = controller1.Ak;
    Bk1 = controller1.Bk;
    Ck1 = controller1.Ck;
    Dk1 = controller1.Dk;

end

%Controller 2

if isempty(controller2)
    
    Ak2=0; Bk2=0; Ck2=0; Dk2=0;
    
else
    
    Ak2 = controller2.Ak;
    Bk2 = controller2.Bk;
    Ck2 = controller2.Ck;
    Dk2 = controller2.Dk;

end

%% SCENARIO 0 - COUPLE OPEN LOOP SYSTEMS FOR OPEN LOOP PERFORMANCE VERIFICATION

%The first scenario couples the two substructures through the
%interconnection stiffness matrix and defines the inputs / outputs as the
%mappings between disturbance inputs and the "virtual" performance outputs

if strcmp(scenario,'scenario_0')==1
    
    %INPUTS:[w1;w2]
    %OUTPUTS:[z1;z2]
    %STATES:[x1;x2]
    

    A = [A1+Bp1*K11*Cq1, Bp1*K12*Cq2;...
        Bp2*K21*Cq1, A2+Bp2*K22*Cq2];
    
    %In the most recent generalized plant formulation, within the function
    %'generalizedplant', the disturbance input, performance output, and of
    %course, the D11 matrices are augmented with additional terms.
    %Tersely, these are related to preventing singularities at high
    %frequencies during controller synthesis.  For examining the open-loop
    %system response, we simply want to look at the [w1;w2] --> [z1;z2]
    %mapping, and so the control input and sensor noise terms that
    %augmented the generalized plant matrices are removed.
    
    Bw1 = Bw1(:,size(Bw1,2)-size(Cy1,1)); %remove additional terms
    Bw2 = Bw2(:,size(Bw2,2)-size(Cy2,1)); %remove additional terms
    
    B = [Bw1, zeros(size(Bw1,1),size(Bw2,2));...
        zeros(size(Bw2,1),size(Bw1,2)), Bw2];
    
    Cz1 = Cz1(size(Cz1,1)-size(Bu1,2),:); %remove additional terms
    Cz2 = Cz2(size(Cz2,1)-size(Bu2,2),:); %remove additional terms
    
    C = [Cz1, zeros(size(Cz1,1),size(Cz2,2));...
        zeros(size(Cz2,1),size(Cz1,2)), Cz2];
    
    D111 = D111(size(Cz1,1),size(D111,2)-size(Cy1,1)); %remove additional terms
    D112 = D112(size(Cz2,1),size(D112,2)-size(Cy2,1)); %remove additional terms
    
    D = [D111, zeros(size(D111,1),size(D112,2));...
        zeros(size(D112,1),size(D111,2)), D112];
    
    %We will now create the system structure.  For this scenario, we are
    %formulating the mapping between all disturbance inputs and performance
    %outputs.
    system = struct('A',A,'B',B,'C',C,'D',D,'ninputs',[],'moutputs',[]);
    
%     disp('INPUTS:[w1;w2], OUTPUTS:[z1;z2], STATES:[x1;x2]')
%     disp('Coupled system has the representation:')
%     disp('[A1+Bp1*K11*Cq1, Bp1*K12*Cq2,     Bw1,    0;')
%     disp('Bp2*K21*Cq1,   A2+Bp2*K22*Cq2,     0,    Bw2;')
%     disp('  Cz1,                 0,         D111,   0')
%     disp('   0,                 Cz2,         0,   D112]')
    
    
end


%% SCENARIO 1 - NO CONTROLLERS, SYNTHESIZE K1

%This scenario poses the problem for controller synthesis around the first
%subsystem.  In this situation, there is no controller for subsystem #2.

if strcmp(scenario,'scenario_1')==1
    
%     disp('Formulating system for controller synthesis --- SCENARIO 1')
    
    %INPUTS:[w1;w2;u1]
    %OUTPUTS:[z1;y1]
    %STATES:[x1;x2]
    
    A = [A1+Bp1*K11*Cq1, Bp1*K12*Cq2;...
        Bp2*K21*Cq1, A2+Bp2*K22*Cq2];
    
    B = [Bw1, zeros(size(Bw1,1),size(Bw2,2)), Bu1;...
        zeros(size(Bw2,1),size(Bw1,2)), Bw2, zeros(size(Bw2,1),size(Bu1,2))];
    
    C = [Cz1 zeros(size(Cz1,1),size(A2,2));...
        Cy1, zeros(size(Cy1,1),size(A2,2))];
    
    D = [zeros(size(Cz1,1),size(Bw1,2)) zeros(size(Cz1,1),size(Bw2,2)) D121;...
        D211 zeros(size(D211,1),size(Bw2,2)) zeros(size(D211,1),size(D121,2))];
    
    %This system structure, which looks like G ~ (A,B,C,D), organizes all
    %inputs and outputs in the following way: 
    %INPUTS = [w1;w2;u1];
    %OUTPUTS = [z1;y1];
    %STATES = [x1;x2];
    %When we use the D-K iteration code, we must input the number of control inputs
    %and outputs to the algorithm.  They are assumed to be the last n
    %inputs and m outputs
    
    n = size(Bu1,2);
    m = size(Cy1,1);
    system = struct('A',A,'B',B,'C',C,'D',D,'ninputs',n,'moutputs',m);
    
%     disp('INPUTS:[w1;w2;u1], OUTPUTS:[z1;y1], STATES:[x1;x2]')
%     disp('Coupled system has the representation:')
%     disp('[A1+Bp1*K11*Cq1, Bp1*K12*Cq2,     Bw1,    0,    Bu1;')
%     disp('Bp2*K21*Cq1,   A2+Bp2*K22*Cq2,     0,    Bw2,    0;')
%     disp('  Cz1,                 0,          0,     0,    D121')
%     disp('  Cy1,                 0,         D211,   0,     0]')
    
end

%% SCENARIO 2 - HAVE K1, SYNTHESIZE K2

if strcmp(scenario,'scenario_2')==1
    
    %INPUTS:[w1;w2,u2]
    %OUTPUTS:[z2;y2]
    %STATES:[x1;x2;xk1]
    
%     disp('Formulating system for controller synthesis --- SCENARIO 2')
    
    A = [A1+Bu1*Dk1*Cy1+Bp1*K11*Cq1, Bp1*K12*Cq2, Bu1*Ck1;...
        Bp2*K21*Cq1, A2+Bp2*K22*Cq2, zeros(size(A2,1),size(Ck1,2));...
        Bk1*Cy1, zeros(size(Bk1,1),size(A2,1)) Ak1];
    
    B = [Bw1+Bu1*Dk1*D211, zeros(size(Bw1,1),size(Bw2,2)) zeros(size(Bw1,1),size(Bu2,2));...
        zeros(size(Bw2,1),size(Bw1,2)), Bw2, Bu2;...
        Bk1*D211, zeros(size(Bk1,1),size(Bw2,2)), zeros(size(Bk1,1),size(Bu2,2))];
    
    C = [zeros(size(Cz2,1),size(A1,2)), Cz2, zeros(size(Cz2,1),size(Ak1,2));...
        zeros(size(Cy2,1),size(A1,2)), Cy2, zeros(size(Cy2,1),size(Ak1,2))];
    
    D = [zeros(size(Cz2,1),size(Bw1,2)), zeros(size(Cz2,1),size(Bw2,2)), D122;...
        zeros(size(Cy2,1),size(Bw1,2)), D212, zeros(size(Cy2,1),size(Bu2,2))];
    
    n = size(Bu2,2);
    m = size(Cy2,1);
    system = struct('A',A,'B',B,'C',C,'D',D,'ninputs',n,'moutputs',m);
   
%     disp('INPUTS:[w1;w2,u2], OUTPUTS:[z2;y2], STATES:[x1;x2;xk1]')
%     disp('Coupled system has the representation:')
%     disp('[A1+Bp1*Dk1*Cy1+Bq1*K11*Cq1, Bp1*K12*Cq2, Bu1*Ck1,   Bw1+Bu1*Dk1*D211,   0,    0;')
%     disp(' Bp2*K21*Cq1,              A2+Bp2*K22*Cq2,   0,             0,          Bw2,  Bu2;')
%     disp('   Bk1*Cu1,                      0,          Ak1,       Bk1*D211,        0,    0;')
%     disp('   0,                           Cz2,         0,             0,           0,   D122')
%     disp('   0,                           Cy2,         0,             0,         D212,   0]')
end

%% SCENARIO 3 - HAVE K2, SYNTHESIZE K1

if strcmp(scenario,'scenario_3')==1
    
    %INPUTS: [w1;w2;u1]
    %OUTPUTS: [z1;y1]
    %STATES: [x1;x2;xk2]
    
%     disp('Formulating system for controller synthesis --- SCENARIO 3')
    
    A = [A1+Bp1*K11*Cq1, Bp1*K12*Cq2, zeros(size(A1,1),size(Ck2,2));...
        Bp2*K21*Cq1, A2+Bp2*K22*Cq2+Bu2*Dk2*Cy2, Bu2*Ck2;...
        zeros(size(Ak2,1),size(A1,2)), Bk2*Cy2, Ak2];
    
    B = [Bw1, zeros(size(Bw1,1),size(Bw2,2)), Bu1;...
        zeros(size(Bw2,1),size(Bw1,2)), Bw2+Bu2*Dk2*D212, zeros(size(Bw2,1),size(Bu1,2));...
        zeros(size(Bk2,1),size(Bw1,2)), Bk2*D212, zeros(size(Bk2,1),size(Bu1,2))];
    
    C = [Cz1 zeros(size(Cz1,1),size(A2,2)), zeros(size(Cz1,1),size(Ak2,2));...
        Cy1, zeros(size(Cy1,1),size(A2,2)), zeros(size(Cy1,1),size(Ak2,2))];
    
    D = [zeros(size(Cz1,1),size(Bw1,2)), zeros(size(Cz1,1),size(Bw2,2)), D121;...
        D211, zeros(size(D211,1),size(Bw2,2)), zeros(size(D211,1),size(D121,2))];
    
    n = size(Bu1,2);
    m = size(Cy1,1);
    system = struct('A',A,'B',B,'C',C,'D',D,'ninputs',n,'moutputs',m);
    
%     disp('INPUTS:[w1;w2,u1], OUTPUTS:[z1;y1], STATES:[x1;x2;xk2]')
%     disp('Coupled system has the representation:')
%     disp('[A1+Bq1*K11*Cp1,         Bq1*K12*Cp2,             0,       Bw1,        0,         Bu1;')
%     disp(' Bq2*K21*Cp1,    A2+Bq2*K22*Cp2+Bu2*Dk2*Cy2,   Bu2*Ck2,     0,  Bw2+Bu2*Dk2*D212,  0;')
%     disp('     0,                      Bk2*Cy2,            Ak2,       0,     Bk2*D212,       0;')
%     disp('   Cz1,                        0,                 0,        0,          0,        D121;')
%     disp('   Cy1,                        0,                 0,       D211,        0,         0]')
    
end

%% SCENARIO 4 - HAVE K1, K2, FORMULATE COMPOSITE SYSTEM FOR PERFORMANCE VERIFICATION

if strcmp(scenario,'scenario_4')==1
    
    %INPUTS: [w1;w2]
    %OUTPUTS: [z1;z2]
    %STATES: [x1;x2;xk1;xk2]
    
%     disp('Formulating system for performance verification --- SCENARIO 4')
    

    A = [A1+Bu1*Dk1*Cy1+Bp1*K11*Cq1, Bp1*K12*Cq2, Bu1*Ck1, zeros(size(A1,1),size(Ak2,2));...
        Bp2*K21*Cq1, A2+Bu2*Dk2*Cy2+Bp2*K22*Cq2, zeros(size(A2,1),size(Ck1,2)), Bu2*Ck2;...
        Bk1*Cy1, zeros(size(Bk1,1),size(A2,2)), Ak1, zeros(size(Ak1,1),size(Ak2,2));...
        zeros(size(Ak2,1),size(A1,2)), Bk2*Cy2, zeros(size(Bk2,1),size(Ak1,2)), Ak2];
    
    B = [Bw1+Bu1*Dk1*D211, zeros(size(Bw1,1),size(Bw2,2));...
        zeros(size(Bw2,1),size(Bw1,2)), Bw2+Bu2*Dk2*D212;...
        Bk1*D211, zeros(size(Bk1,1),size(Bw2,2));...
        zeros(size(Bk2,1),size(Bw1,2)), Bk2*D212];
    
    C = [Cz1+D121*Dk1*Cy1, zeros(size(Cz1,1),size(Cz2,2)), D121*Ck1, zeros(size(Cz1,1),size(Ck2,2));...
        zeros(size(Cz2,1),size(Cz1,2)), Cz2+D122*Dk2*Cy2, zeros(size(Cz2,1),size(Ck1,2)), D122*Ck2];
    
    D = [D121*Dk1*D211, zeros(size(D121,1),size(D212,2));...
        zeros(size(D122,1),size(D211,2)), D122*Dk2*D212];
   
    n = [];
    m = [];
    system = struct('A',A,'B',B,'C',C,'D',D,'ninputs',n,'moutputs',m);
    
%     disp('INPUTS:[w1;w2], OUTPUTS:[z1;y2], STATES:[x1;x2;xk1;xk2]')
%     disp('Coupled system has the representation:')
%     disp('[A1+Bu1*Dk1*Cy1+Bp1*K11*Cq1,     Bp1*K12*Cq2,        Bu1*Ck1,     0,    Bw1+Bu1*Dk1*D211,          0;')
%     disp('     Bp2*K21*Cq1,       A2+Bu2*Dk2*Cy2+Bp2*K22*Cq2,     0,     Bu2*Ck2,        0,          Bw2+Bu2*Dk2*D212;')
%     disp('       Bk1*Cy1,                      0,                Ak1,       0,        Bk1*D211,              0;')
%     disp('          0,                      Bk2*Cy2,              0,       Ak2,          0,               Bk2*D212;')
%     disp('   Cz1+D121*Dk1*Cy1,                 0,             D121*Ck1,     0,      D121*Dk1*D211,           0;')
%     disp('          0,                 Cz2+D122*Dk2*Cy2,          0,     D122*Ck2,        0,          D122*Dk2*D212]')  
    
end

if strcmp(scenario,'scenario_5')==1
    
    %INPUTS: [w1;w2]
    %OUTPUTS: [z1;z2]
    %STATES: [x1;x2;xk1;xk2]
    
    %     disp('Formulating system for performance verification with D12 and D21 terms zeroed --- SCENARIO 5')

    D121 = 0*D121;
%     D211 = 0*D211;
    D122 = 0*D122;
%     D212 = 0*D212;
    
    A = [A1+Bu1*Dk1*Cy1+Bp1*K11*Cq1, Bp1*K12*Cq2, Bu1*Ck1, zeros(size(A1,1),size(Ak2,2));...
        Bp2*K21*Cq1, A2+Bu2*Dk2*Cy2+Bp2*K22*Cq2, zeros(size(A2,1),size(Ck1,2)), Bu2*Ck2;...
        Bk1*Cy1, zeros(size(Bk1,1),size(A2,2)), Ak1, zeros(size(Ak1,1),size(Ak2,2));...
        zeros(size(Ak2,1),size(A1,2)), Bk2*Cy2, zeros(size(Bk2,1),size(Ak1,2)), Ak2];
    
    B = [Bw1+Bu1*Dk1*D211, zeros(size(Bw1,1),size(Bw2,2));...
        zeros(size(Bw2,1),size(Bw1,2)), Bw2+Bu2*Dk2*D212;...
        Bk1*D211, zeros(size(Bk1,1),size(Bw2,2));...
        zeros(size(Bk2,1),size(Bw1,2)), Bk2*D212];
    
    C = [Cz1+D121*Dk1*Cy1, zeros(size(Cz1,1),size(Cz2,2)), D121*Ck1, zeros(size(Cz1,1),size(Ck2,2));...
        zeros(size(Cz2,1),size(Cz1,2)), Cz2+D122*Dk2*Cy2, zeros(size(Cz2,1),size(Ck1,2)), D122*Ck2];
    
    D = [D121*Dk1*D211, zeros(size(D121,1),size(D212,2));...
        zeros(size(D122,1),size(D211,2)), D122*Dk2*D212];
   
    n = [];
    m = [];
    
    system = struct('A',A,'B',B,'C',C,'D',D,'ninputs',n,'moutputs',m);

    
end


end