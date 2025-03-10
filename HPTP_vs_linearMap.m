clear;

%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com

v0 = [1;0;0];
v1 = [0;1;0];
v2 = [0;0;1];
G3 = 1/sqrt(3)*(Tensor(v0,v0,v0) + Tensor(v1,v1,v1) + Tensor(v2,v2,v2));
da = 3; db=3; dc=3;

fHP = [];
fHPTP = [];

for p=0:0.1:1
rhoABC = (1-p)*(G3*G3') + p*eye(da^3)/da^3; % initial state
rhoAB = PartialTrace(rhoABC, 3, [da db dc]);

% recovery by linear maps
cvx_begin sdp quiet
    cvx_solver sedumi
    variable J(db*db*dc, db*db*dc) complex
    
    sigABC = PartialTrace(kron(eye(da), J) * kron(PartialTranspose(rhoAB, 2, [da db]), eye(db*dc)), 2, [da db db*dc]);
    f = TraceNorm(sigABC-rhoABC);
    minimize f
cvx_end
fHP(end+1) = f;

% recovery by HPTP maps
cvx_begin sdp quiet
    cvx_solver sedumi
    variable J(db*db*dc, db*db*dc) hermitian
   
    sigABC = PartialTrace(kron(eye(da), J) * kron(PartialTranspose(rhoAB, 2, [da db]), eye(db*dc)), 2, [da db db*dc]);
    f = TraceNorm(sigABC-rhoABC);
    minimize f 
    subject to
        PartialTrace(J, 2, [db db*dc]) == eye(db);
cvx_end
fHPTP(end+1) = f;
end