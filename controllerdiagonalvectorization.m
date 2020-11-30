

function [Avecs,Bvecs,Cvecs,Dvecs,nk1,mk1,rk1,numrealeigs,index_k] = controllerdiagonalvectorization(Controller,tol)


%Unpackage controller
K1 = Controller;
Ak1 = K1.A;
Bk1 = K1.B;
Ck1 = K1.C;
Dk1 = K1.D;

%Get some matrix sizes
nk1 = size(Ak1,1);
mk1 = size(Bk1,2);
rk1 = size(Ck1,1);

%Find controller eigenvectors
[Tk1,~] = eig(Ak1);

%Convert controller into modal form
Ak1m = Tk1\Ak1*Tk1;
Bk1m = Tk1\Bk1;
Ck1m = Ck1*Tk1;
Dk1m = Dk1;

Ak1m_diag = diag(Ak1m);

%Identify the truly real eigenvalues

% tol = 1e-3;

Ak1m_diagreal = real(Ak1m_diag);
Ak1m_diagimag = imag(Ak1m_diag);

index_k1 = [];
for i = 1:length(Ak1m_diag)

zeta = abs(Ak1m_diagreal(i)/abs(Ak1m_diag(i)));
    if zeta >= 1-tol
        Ak1m_diagimag(i) = 0;
        index_k1 = [index_k1 i];
    end
    
end

%This means that the eigenvector corresponding to these eigenvalues is also
%real.  So we can remove the imaginary components in the B and C matrices
%that correspond to this.  This means we are looking to the row index in B
%and the columns in C.

Bk1m_real = real(Bk1m);
Bk1m_imag = imag(Bk1m);
Ck1m_real = real(Ck1m);
Ck1m_imag = imag(Ck1m);

Bk1m_imag(index_k1,:) = 0;
Ck1m_imag(:,index_k1) = 0;


%% State reordering
%Ultimately, we need a way to track the real-valued eigenvalues for
%reconstruction.  Global variable creation between main script and fitness
%function may be a way to achieve this end.

%A matrix
Ak1m_realzero = Ak1m_diagreal(index_k1,:);
Ak1m_imagzero = Ak1m_diagimag(index_k1,:);

Ak1m_realzerosremoved = Ak1m_diagreal;
Ak1m_imagzerosremoved = Ak1m_diagimag;

Ak1m_realzerosremoved(index_k1,:) = [];
Ak1m_imagzerosremoved(index_k1,:) = [];

Ak1m_realreordered = [Ak1m_realzero;...
    Ak1m_realzerosremoved];
Ak1m_imagreordered = [Ak1m_imagzero;...
    Ak1m_imagzerosremoved];

%B matrix
Bk1m_realzero = Bk1m_real(index_k1,:);
Bk1m_imagzero = Bk1m_imag(index_k1,:);

Bk1m_realzerosremoved = Bk1m_real;
Bk1m_imagzerosremoved = Bk1m_imag;

Bk1m_realzerosremoved(index_k1,:) = [];
Bk1m_imagzerosremoved(index_k1,:) = [];

Bk1m_realreordered = [Bk1m_realzero; Bk1m_realzerosremoved];
Bk1m_imagreordered = [Bk1m_imagzero; Bk1m_imagzerosremoved];

%C matrix
Ck1m_realzero = Ck1m_real(:,index_k1);
Ck1m_imagzero = Ck1m_imag(:,index_k1);

Ck1m_realzerosremoved = Ck1m_real;
Ck1m_imagzerosremoved = Ck1m_imag;

Ck1m_realzerosremoved(:,index_k1) = [];
Ck1m_imagzerosremoved(:,index_k1) = [];

Ck1m_realreordered = [Ck1m_realzero, Ck1m_realzerosremoved];
Ck1m_imagreordered = [Ck1m_imagzero, Ck1m_imagzerosremoved];

numrealeigs = length(index_k1);

Avecs = [Ak1m_realreordered, Ak1m_imagreordered];
Bvecs = [Bk1m_realreordered, Bk1m_imagreordered];
Cvecs = [Ck1m_realreordered;Ck1m_imagreordered];
Dvecs = Dk1;

index_k = index_k1;

end