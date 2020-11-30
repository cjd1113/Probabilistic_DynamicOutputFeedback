

function [Mcc,Kcc,KE,ME,Le,PD,PC,Mm] = bernoullibeamFEMfuncboundarydef(E,rho,b,h,Ne,TotalLength,Nc,Nd,Meas,constraintloc)

%Finite element modeling of an Euler Bernoulli Beam
%Prepared by: Christopher D'Angelo
%Date: October 26, 2017

%Generate beam model 
%final value of 1 means beginning of beam is constrained, 2 means the end
%is constrained, 3 means pinned-pinned, and 0 means unconstrained

A = b*h;
I = 1/12 * b*h^3;

L = TotalLength/Ne;

Ke = (E*I)/(L^3) * [12 6*L -12 6*L;...
    6*L 4*L^2 -6*L 2*L^2;...
    -12 -6*L 12 -6*L;...
    6*L 2*L^2 -6*L 4*L^2];

Me = (rho*A*L/420)* [156 22*L 54 -13*L;...
    22*L 4*L^2 13*L -3*L^2;...
    54 13*L 156 -22*L;...
    -13*L -3*L^2 -22*L 4*L^2];

Nodes = Ne*2+2; % Number of Nodes in the model - rotational, transverse

%Assembly via the direct stiffness method
Le = zeros(4,Nodes,Ne); %Assembly matrices
I = eye(4); %Each element has 4 dof

%%% Mass and Stiffness element matrices, along with element locator %%%

KE = zeros(Nodes,Nodes,Ne);
ME = zeros(Nodes,Nodes,Ne);

for i = 1:Ne
    if i==1
        Le(:,:,i) = [I zeros(4,Nodes-4)];
        KE(:,:,i) = Le(:,:,i)'*Ke*Le(:,:,i);
        ME(:,:,i) = Le(:,:,i)'*Me*Le(:,:,i);
    end
    if i==2
        Le(:,:,i) = [zeros(4,2) I zeros(4,Nodes-4-2)];
        KE(:,:,i) = Le(:,:,i)'*Ke*Le(:,:,i);
        ME(:,:,i) = Le(:,:,i)'*Me*Le(:,:,i);
        inc = 2;
    end
    if i>2
        inc = inc+2;
        Le(:,:,i) = [zeros(4,inc) I zeros(4,Nodes-4-inc)];
        KE(:,:,i) = Le(:,:,i)'*Ke*Le(:,:,i);
        ME(:,:,i) = Le(:,:,i)'*Me*Le(:,:,i);
        
    end
end

KK = sum(KE,3);
MM = sum(ME,3);

%%% DISTURBANCE AND CONTROL INPUT LOCATIONS %%%

%%Now onto load vectors.  

Pe = [L/2; 1/12*L^2; L/2; -1/12*L^2];
f = 1; %without loss of generality, force is equal to 1
Pe = f*Pe;

%Separate it into disturbances and control inputs

%Disburbance elements

if any(Nd>Ne)==1
    disp('Disturbance element locations fall outside feasible range')
    return
end
if any(Nc>Ne)==1
    disp('Control element locations fall outside feasible range')
    return
end

% Nd = [3 14]; %Disturbances enter elements 3 and 14
% %Control input elements
% Nc = [6 18]; %Control inputs enter elements 6 and 18

% Pd = zeros(Nodes,1,Ne);
% Pc = zeros(Nodes,1,Ne);
% 
% for i = 1:Ne
%     if any(i==Nd)==1
%         Pd(:,:,i) = Le(:,:,i)'*Pe;
%     end
%     
%     if any(i==Nc)==1
%         Pc(:,:,i) = Le(:,:,i)'*Pe;
%     end
% end
% 
% PD = sum(Pd,3);
% PC = sum(Pc,3);

dinputs = length(Nd);
cinputs = length(Nc);

Pd = zeros(Nodes,1,Ne,dinputs);
Pc = zeros(Nodes,1,Ne,cinputs);

for i = 1:Ne
    for j = 1:length(Nd)
        if Nd(j)==i
            Pd(:,:,i,j) = Le(:,:,i)'*Pe;
        end
        if Nc(j)==i
            Pc(:,:,i,j) = Le(:,:,i)'*Pe;
        end
    end
    
end

PD = sum(Pd,3);
PC = sum(Pc,3);

for i = 2:size(PD,4)
   PD = [PD(:,:,:,i-1) PD(:,:,:,i)];
end

for i = 2:size(PC,4)
    PC = [PC(:,:,:,i-1) PC(:,:,:,i)];
end


%%% MEASUREMENT LOCATIONS %%%

%%Odd terms are transverse displacements, even terms are rotations

%%Let us define some measurements, only in the transverse direction

Me = [1;0;1;0];
% Nm = [5 50 75]; %measurements made at elements 5, 50, 75
Nm = Meas;
if any(Nm>Ne)==1
    disp('Measurement locations fall outside of feasible range')
    return
end
% Mm = zeros(length(Nm),Nodes,Ne);
Mm = zeros(1,Nodes);
for i = 1:Ne
    if any(i==Nm)==1
        Mmtemp = Me'*Le(:,:,i);
        Mm = [Mm;Mmtemp];
    end
end
Mm = Mm(2:end,:);


%%% CONSTRAINTS %%%
%Now, let's constrain one end

if constraintloc == 1

    Mcc = MM(3:end,3:end);
    Kcc = KK(3:end,3:end);
    PD = PD(3:end,:);
    PC = PC(3:end,:);
    Mm = Mm(:,3:end);
    
end

if constraintloc == 2

    Mcc = MM(1:end-2,1:end-2);
    Kcc = KK(1:end-2,1:end-2);
    PD = PD(1:end-2,:);
    PC = PC(1:end-2,:);
    Mm = Mm(:,1:end-2);
    
end

if constraintloc == 3

    Mcc = MM(3:end-2,3:end-2);
    Kcc = KK(3:end-2,3:end-2);
    PD = PD(3:end-2,:);
    PC = PC(3:end-2,:);
    Mm = Mm(:,3:end-2);
    
end

if constraintloc == 0
    
    Mcc = MM;
    Kcc = KK;
    PD = PD;
    PC = PC;
    Mm = Mm;
    
end

return
