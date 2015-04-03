 % Computing the quivalent MIL repsentation of STL formulae

function [Matrices,z,lbz,ubz] = get_MILP_matrices(A,B,C,nx,nu,nw,Wref,wlb,wub,Np,M)

% %testing the file
% A = rand(5,5);
% B = rand(5,4);
% C = rand(5,7);
% Wref = rand(7,1);
% wub = 1;
% wlb = -1;
% Np = 3;
% 
% nx = size(A,2);
% nu = size(B,2);
% nw = size(C,2);
nz = 3;
% M = 1000;
% state-space model x(t) = A x(t-1) + B u(t-1) + C w(t-1)
% E1*x(k) + E2*u(k) + E3*z(k) + E4*w(k) <= g5
% H*z(k) = f

E1 = [zeros(1,nx);
      A(5,:);
      zeros(2,nx);
      A(5,:);
      -A(5,:)
      ];
E2 = [zeros(1,nu);
      B(5,:);
      zeros(2,nu);
      B(5,:);
      -B(5,:)
      ];
E3 = [-1 0 0;
      -1 0 0;
      -1 M 0;
      1 M 0;
      -1 0 M;
      1 0 M
      ];
E4 = [zeros(1,nw-1) -1000;
      C(5,:)+[zeros(1,5) -1 0];
      [zeros(2,nw-1) [-1000;1000]];
      C(5,:)+[zeros(1,5) -1 0];
      -C(5,:)+[zeros(1,5) 1 0]
      ];
g5 = [0;0;M;M;M;M];
H = [0 1 1];
f = 1;

L=size(E1,1);           % number of inequalities

S = [1 zeros(1,nw-1);
    0 1 zeros(1,nw-2); 
    0 0 1 zeros(1,nw-3),
    zeros(1,3) 1 0 0 0;
    zeros(1,4) 1 0 0;
    zeros(1,5) 1 0;
    zeros(1,6) 1;
    -1 zeros(1,nw-1);
    0 -1 zeros(1,nw-2);
    0 0 -1 zeros(1,nw-3);
    zeros(1,3) -1 0 0 0;
    zeros(1,4) -1 0 0;
    zeros(1,5) -1 0;
    zeros(1,6) -1
    ];               % matrix S in S w <= q
q1 = [Wref;-Wref]+ [wub'*ones(1,size(Wref,2));-wlb'*ones(1,size(Wref,2))];
q = q1';            % to obtain the constraint vector having each horizon on one row
% defining the matrices for MPC problem
% Cmpc x(k) + Dmpc \tilde v(k) + Rmpc \tilde w(k) <= Empc
% Hmpc \tilde z(k) = Fmpc 
Cmpc=[];
Empc = [];
Hmpc = [];

for i=1:2*Np
   Cmpc = [Cmpc; E1*A^i];
   Empc = [Empc; g5];
   Hmpc = [Hmpc H];
end

p1 = E2;
p2 = E3; 
p3 = E4;
for j=0:2*Np-2
    p1=[p1; E1*A^j*B];
    p3=[p3; E1*A^j*C];
end
p2 = [p2;zeros(size(p1,1)-L,size(E3,2))];
ii1 = 0;
ii2 = 0;
ii3 = 0;
for i=0:2*Np-1
    Dmpc1(:,ii1+1:ii1+nu)=[zeros(i*L,nu); p1(1:L*(2*Np-i),:)];
    Dmpc2(:,ii2+1:ii2+nz)=[zeros(i*L,nz); p2(1:L*(2*Np-i),:)];
    Dmpc3(:,ii3+1:ii3+nw)=[zeros(i*L,nw); p3(1:L*(2*Np-i),:)];
    ii1 = ii1+nu;
    ii2 = ii2 + nz;
    ii3 = ii3 + nw;
end
Dmpc = [Dmpc1 Dmpc2];
Rmpc = Dmpc3;
Hmpc = [zeros(1,nu*(2*Np-1)) Hmpc];
% Fmpc = [zeros(nu*(2*Np-1),1);f*ones(nz*2*Np,1)];
Fmpc = f;

% building the S matrixx for the entire prediction horizon
Sp = [S;zeros(2*nw*(2*Np-1),nw)]; 
iie = 0;
Ls = size(S,1);
for ie=0:2*Np-1
     Smpc(:,iie+1:iie+nw)=[zeros(ie*Ls,nw); Sp(1:Ls*(2*Np-ie),:)];
    iie = iie+nw;
end
qmpc = reshape(q.',[],1);
%% Objective function
% J(k) = C1(k)^T(k)u(k) + C2(k)^T(k) z(k) + C4^T(k)w(k)

C1 = zeros(nu*2*Np,1);
C2 = [1; zeros(nz*2*Np-1,1)];
C3 = zeros(nw*2*Np,1);
CC = -[C1;C2];
%% bounds on z (STL) variables
z = [sdpvar(1,1);binvar(2,1)];
lbz = [0.000001;0;0];
ubz = [Inf;1;1];
for kz=1:2*Np-1
   z = [z; [sdpvar(1,1);binvar(2,1)]]; 
   lbz = [lbz; [0.000001;0;0]];
   ubz = [ubz; [1000;1;1]];
end

%% defining the output matrices

Matrices.Cmpc = Cmpc;
Matrices.Dmpc1 = Dmpc1;
Matrices.Dmpc2 = Dmpc2;
Matrices.Rmpc = Rmpc;
Matrices.Empc = Empc;
Matrices.Hmpc = Hmpc;
Matrices.Fmpc = Fmpc;
Matrices.Smpc = Smpc;
Matrices.qmpc = qmpc;
Matrices.CC1 = C1;
Matrices.CC2 = -C2;
