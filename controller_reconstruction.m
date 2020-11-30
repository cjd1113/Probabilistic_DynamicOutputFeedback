



function [Ak,Bk,Ck,Dk] = ...
    controller_reconstruction(ControllerVariables,numrealeigs,nk,mk,rk)


%Prepared by: Chris D'Angelo
%Date: August 27, 2018

%Remember that the zero imaginary parts were put at the beginning of each
%imaginary vector.

%Controller order is nk.  This means we have nk eigenvalues.  We could have
%nk real-valued eigenvalues, or nk complex valued eigenvalues.  If we have
%nk complex valued eigenvalues, we have nk/2 real components and nk/2
%complex valued components.  Thus, the first (nk/2 + numrealeigs/2 elements
%will correspond to the real-parts of the eigenspectrum.

%Check if row vector, if true, make column.  This is to match how Matlab's
%ga function handles decision variables --- as column vectors.
if isrow(ControllerVariables)==1
    ControllerVariables = ControllerVariables';
end

%%START WITH CONSTRUCTION OF AK
Akreal = ControllerVariables(1:(nk/2+numrealeigs/2)); %separate out the real parts
Akimag = [zeros(numrealeigs,1);...
    ControllerVariables(length(Akreal)+1:length(Akreal)+nk/2-numrealeigs/2)];

Akrealimag = [Akreal Akimag];

%%NOW ONTO CONSTRUCTION OF BK
%We are now at nk+1 along the ControllerVariables vector.  Staring now,
%with Bk construction

Bkreal = ControllerVariables(nk+1:nk+mk*(nk/2+numrealeigs/2));
Bkrealmatrix = reshape(Bkreal,nk/2+numrealeigs/2,mk);

Bkimag = ControllerVariables(nk+mk*(nk/2+numrealeigs/2)+1:...
    nk+mk*(nk/2+numrealeigs/2)+(mk*(nk/2-numrealeigs/2)));
Bkimagmatrixtemp = reshape(Bkimag,nk/2-numrealeigs/2,mk);
Bkimagmatrix = [zeros(numrealeigs,mk);...
    Bkimagmatrixtemp];

% Bkrealimag = [Bkrealmatrix Bkimagmatrix];

%%NOW ONTO CONSTRUCTION OF CK

Ckbegin = nk+mk*(nk/2+numrealeigs/2)+(mk*(nk/2-numrealeigs/2))+1;

Ckreal = ControllerVariables(Ckbegin:Ckbegin-1+rk*(nk/2+numrealeigs/2));

Ckimag = ControllerVariables(Ckbegin+rk*(nk/2+numrealeigs/2):end);

Ckrealmatrix = reshape(Ckreal,nk/2+numrealeigs/2,rk);
Ckimagmatrixtemp = reshape(Ckimag,nk/2-numrealeigs/2,rk);
Ckimagmatrix = [zeros(numrealeigs,rk);...
    Ckimagmatrixtemp];

% Ckrealimag = [Ckrealmatrix Ckimagmatrix];


%%LEFT OFF HERE - 8/24/18
%% NEXT STEP IS TO TRANSFORM THESE CONTROLLER VARIABLES BACK INTO THEIR COMPLEX STATE SPACE FORM...THEN GO!!!!

%%%%%% Ak %%%%%%

%First step is to find and separate out all of the complex components.
Akrealcomplex = Akrealimag(numrealeigs+1:end,1);
Akcomplexcomplex = Akrealimag(numrealeigs+1:end,2);
%Handle the real parts first
Arealtemp_new = kron(Akrealcomplex',[1 0]);
Arealtemp_new1 = Arealtemp_new';
Arealtemp_new2 = [0;Arealtemp_new(1:end-1)'];
Arealtemp = Arealtemp_new1 + Arealtemp_new2;
%Now onto the complex parts
Acomplextemp_new = kron(Akcomplexcomplex',[1,0]);
Acomplextemp_new1 = Acomplextemp_new';
Acomplextemp_new2 = -[0;Acomplextemp_new(1:end-1)'];
Acomplextemp = Acomplextemp_new1 + Acomplextemp_new2;

Arealcomplex = Arealtemp + 1j.*Acomplextemp;

Arealcomplextotal = [Akrealimag(1:numrealeigs,1) + 1j.*zeros(numrealeigs,1);...
    Arealcomplex];

Ak = eye(nk).*Arealcomplextotal;

%%%%%% Bk %%%%%%

% Bkrealimag = [Bkrealmatrix Bkimagmatrix];

Bkrealcomplex = Bkrealmatrix(numrealeigs+1:end,:);
Bkcomplexcomplex = Bkimagmatrix(numrealeigs+1:end,:);

Bkrealcomplex_temp1 = kron(Bkrealcomplex',[1 0]);
Bkrealcomplex_temp2 = [zeros(mk,1), Bkrealcomplex_temp1(:,1:end-1)];
Bkrealcomplex = Bkrealcomplex_temp1' + Bkrealcomplex_temp2';

Bkcomplexcomplex_temp1 = kron(Bkcomplexcomplex',[1 0]);
Bkcomplexcomplex_temp2 = -[zeros(mk,1),Bkcomplexcomplex_temp1(:,1:end-1)];
Bkcomplexcomplex = Bkcomplexcomplex_temp1' + Bkcomplexcomplex_temp2';

Bkrealcomplex = Bkrealcomplex + 1j.*Bkcomplexcomplex;

Bkrealcomplextotal = [Bkrealmatrix(1:numrealeigs,:) + 1j.*zeros(numrealeigs,mk);...
    Bkrealcomplex];

Bk = Bkrealcomplextotal;

%%%%%% Ck %%%%%%

% Ckrealimag = [Ckrealmatrix Ckimagmatrix]

Ckrealcomplex = Ckrealmatrix(numrealeigs+1:end,:);
Ckcomplexcomplex = Ckimagmatrix(numrealeigs+1:end,:);

Ckrealcomplex_temp1 = kron(Ckrealcomplex',[1 0]);
Ckrealcomplex_temp2 = [zeros(rk,1),Ckrealcomplex_temp1(:,1:end-1)];
Ckrealcomplex = Ckrealcomplex_temp1' + Ckrealcomplex_temp2';

Ckcomplexcomplex_temp1 = kron(Ckcomplexcomplex',[1 0]);
Ckcomplexcomplex_temp2 = -[zeros(rk,1), Ckcomplexcomplex_temp1(:,1:end-1)];
Ckcomplexcomplex = Ckcomplexcomplex_temp1' + Ckcomplexcomplex_temp2';

Ckrealcomplex = Ckrealcomplex + 1j.*Ckcomplexcomplex;

Ckrealcomplextotal = [Ckrealmatrix(1:numrealeigs,:) + 1j.*zeros(numrealeigs,rk);...
    Ckrealcomplex];

Ck = Ckrealcomplextotal.';

%%%%%% Dk %%%%%%

Dk = zeros(rk,mk);

% ControllerVariables - (Ak,Bk,Ck,Dk)


end