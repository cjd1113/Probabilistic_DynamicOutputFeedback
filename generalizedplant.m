
function [A, B1, B2, C1, C2, D11, D12, D21] = ...
    generalizedplant(Mcc,Ccc,Kcc,PD,PC,Mm,V,D,filteroption,coordinates)



if strcmp(filteroption,'filter')==1
    
%     [Ccc,V,D] = addDAMPING(Mcc,Kcc,0.01);
    
    if strcmp(coordinates, 'physical')==1
       sprintf('Coordinates must be modal for this option')
       return
    else
       
    
        [Mmodal,Cmodal,Kmodal,Pwmodal,Pumodal,Amodal,...
            B1modal,B2modal,C1modal,C2modal] = ...
            modalcoords(Mcc,Ccc,Kcc,PD,PC,Mm,V,D);
        
        wc = 500;
        [alpha] = filterfunctions(Amodal,B1modal,C1modal,D,wc);
        
        B1modalfilt = [B1modal(1:size(B1modal,1)/2,:);...
            B1modal(size(B1modal,1)/2+1:end,:).*alpha];
        
        [n,m] = size(B2modal);
        [p,q] = size(C1modal);
        [s,t] = size(C2modal);
        
        r12 = 1e-5;
        D12 = r12*[zeros(p,m);...
            eye(m)];
        
        r21 = 1e-5;
        D21 = r21*[zeros(s,size(B1modal,2)), eye(s)];
        
        C1modalaug = [C1modal;...
            zeros(size(D12,1)-p,t)];
        
        D11 = zeros(size(C1modalaug,1),size(D21,2));
        
        B1modalaug = [B1modalfilt, zeros(size(B1modalfilt,1),...
            size(D11,2)-size(B1modalfilt,2))];
        
        A = Amodal;
        B1 = B1modalaug;
        B2 = B2modal;
        C1 = C1modalaug;
        C2 = C2modal;
    
    end
    
elseif strcmp(filteroption,'nofilter')==1
    
%     [Ccc,V,D] = addDAMPING(Mcc,Kcc,0.01);
    
    if strcmp(coordinates,'modal')==1
    
        
        [Mmodal,Cmodal,Kmodal,Pwmodal,Pumodal,Amodal,...
            B1modal,B2modal,C1modal,C2modal] = ...
            modalcoords(Mcc,Ccc,Kcc,PD,PC,Mm,V,D);
        
        alpha = 1;
        
        B1modalfilt = [B1modal(1:size(B1modal,1)/2,:);...
            B1modal(size(B1modal,1)/2+1:end,:).*alpha];
        
        [n,m] = size(B2modal);
        [p,q] = size(C1modal);
        [s,t] = size(C2modal);
        
        r12 = 1e-5;
        D12 = r12*[zeros(p,m);...
            eye(m)];
        
        r21 = 1e-5;
        D21 = r21*[zeros(s,size(B1modal,2)), eye(s)];
        
        C1modalaug = [C1modal;...
            zeros(size(D12,1)-p,t)];
        
        D11 = zeros(size(C1modalaug,1),size(D21,2));
        
        B1modalaug = [B1modalfilt, zeros(size(B1modalfilt,1),...
            size(D11,2)-size(B1modalfilt,2))];
        
        A = Amodal;
        B1 = B1modalaug;
        B2 = B2modal;
        C1 = C1modalaug;
        C2 = C2modal;
    
    elseif strcmp(coordinates,'physical')==1
        
        [Amodal,B1modal,B2modal,C1modal,C2modal,~,~,~] = ...
            physicalcoords(Mcc,Ccc,Kcc,PD,PC,Mm);
        
        alpha = 1;
        
        B1modalfilt = [B1modal(1:size(B1modal,1)/2,:);...
            B1modal(size(B1modal,1)/2+1:end,:).*alpha];
        
        [n,m] = size(B2modal);
        [p,q] = size(C1modal);
        [s,t] = size(C2modal);
        
        r12 = 1e-5;
        D12 = r12*[zeros(p,m);...
            eye(m)];
        
        r21 = 1e-5;
        D21 = r21*[zeros(s,size(B1modal,2)), eye(s)];
        
        C1modalaug = [C1modal;...
            zeros(size(D12,1)-p,t)];
        
        D11 = zeros(size(C1modalaug,1),size(D21,2));
        
        B1modalaug = [B1modalfilt, zeros(size(B1modalfilt,1),...
            size(D11,2)-size(B1modalfilt,2))];
        
        A = Amodal;
        B1 = B1modalaug;
        B2 = B2modal;
        C1 = C1modalaug;
        C2 = C2modal;
        
    end
    
end


end
