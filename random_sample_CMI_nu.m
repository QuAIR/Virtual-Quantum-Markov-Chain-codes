clear;

%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com

%% prepare states
da = 2;
db = 2;
dc = 2;

v0 = [1;0];
v1 = [0;1];

sample_n = 50;


cmi_range = 0.05:0.05:1.6;
len = length(cmi_range);
cmi_list = [];
cost_list = [];

cnt = 1;
for interval = 1:len-1
    i = 1;
while i < sample_n + 1
random_rank = randi([1 8]);
rhoABC = RandomDensityMatrix(da*db*dc, 0, random_rank);  % Random tripartite state

cmi_rho = CMI(rhoABC);
if cmi_rho > cmi_range(interval+1)
    continue
end
if cmi_rho < cmi_range(interval)
    continue
end

rhoAB = PartialTrace(rhoABC, 3, [da db dc]);


%% HPTP recover rho 

JI = MaxEntangled(2,0,1)*MaxEntangled(2,0,1)'*eye(4)*da;  % identity map
cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
    variable JD1(db*db*dc, db*db*dc) hermitian
    variable JD2(db*db*dc, db*db*dc) hermitian
    variable p1
    variable p2
    JD = JD1 - JD2;
    cost = p1 + p2;

    JN = PermuteSystems(kron(JI, JD), [1 3 2 4], [da da db db*dc]);
    sigABC = ApplyMap(rhoAB, JN);
  
    minimize cost
    subject to
        JD1 >= 0; 
        PartialTrace(JD1, 2, [db db*dc]) == p1 * eye(db);
        JD2 >= 0; 
        PartialTrace(JD2, 2, [db db*dc]) == p2 * eye(db);
        sigABC == rhoABC; 

cvx_end


cost_list(cnt) = log2(cost);
cmi_list(cnt) = cmi_rho;
i = i + 1;
cnt = cnt + 1;
end
end
scatter([cmi_list, 0], [cost_list, 0], 'x')
xlabel('CMI','interpreter','latex')
ylabel('$\nu(\rho_{ABC})$', 'interpreter','latex')
% 
save('sample_0.05.mat','cmi_list', 'cost_list')