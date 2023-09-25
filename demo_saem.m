function demo_saem(sim)
addpath(genpath(strcat(pwd,'/MVAR_simu')));
addpath(genpath(strcat(pwd,'/joint_est')));
addpath(strcat(pwd,'/MVAR_simu/data'))
meas_snr=5;
bio_snr_vals=[1 3 5 10];
num_sources=3;
bio_snr_arr=repmat(bio_snr_vals, [50 1]);
dist_thr=30;
T=300;
num_bio_ns=500;
included_sp=1:2004;



%load head model
tmp=load('sim_data');
mesh_=tmp.mesh;
G=tmp.G;
%make results reproducible
[status, seed] = system('od /dev/urandom --read-bytes=4 -tu | awk ''{print $2}''');
seed=str2double(seed);
rng(seed);
sd = rng;

[source_inds] = select_source_inds(mesh_, num_sources, dist_thr);
type=0;
[q, A] = three_source_model(T,type);
[Y, sigma_b, sigma_m, est_snr] = gen_MEG_simu(q,G,source_inds, meas_snr, bio_snr_arr(sim), num_bio_ns);

%SAEM iteration parameters
num_iter = 1000;
b=199;
kappa = 1; 
gamma = zeros(1,num_iter);
gamma(1:2) = 1;
gamma(3:b) = 0.98;
gamma(b+1:end) = 0.98*(((0:num_iter-(b+1))+kappa)/kappa).^(-0.7);

%set number of sources, particles, order etc
Ns=3;
Np=500;
P=2;
scale=1;
whiten_flag=0;
dist_thr = 10; %in mm, min distance between source dipoles in each particle


%initial value for A and V
mean_A0 =zeros(1,Ns)-0.9;
Sigma_A0 = 1e-3*eye(Ns);
lambdamax=10;
while lambdamax >1 || lambdamax < 0.9       
for jj = 1:P
    A_p = mean_A0' + chol(Sigma_A0)*randn(Ns,1);         %diagonal entries of A matrix are U[a,b]
    A_p = diag(diag(A_p));
    A0(Ns*(jj-1)+1:Ns*(jj-1)+Ns,1:Ns) = A_p.*eye(Ns);
end 
A0 = A0';
A0(Ns+1:Ns*P, 1:Ns*(P-1)) = eye(Ns*(P-1));
A0(Ns+1:Ns*P, 1+Ns*(P-1):end) = zeros(Ns*(P-1), Ns);
lambda=eig(squeeze(A0));lambdamax=max(abs(lambda));
end  

V_init_range = [0.1 0.9];
V_q0 = diag(V_init_range(1)+(V_init_range(2)-V_init_range(1))*rand(Ns*P,1));
V_q0=eye(Ns*P)  ;
V_q0(Ns+1:Ns*P, 1:Ns)=0;
V_q0(1:Ns*P, Ns+1:end)=0;

% init struct for SAEM algorithm
init.A0 =A0;
init.V0 = V_q0;
init.q0 = zeros(Ns*P,1);
init.P0=1*eye(Ns*P);
init.sigma_m0 = 1; %eye(size(G,1));
init.sigma_b0 = 5;

% model struct for SAEM
model.G=double(G)*scale;
model.mesh=mesh_;
model.incl_verts=included_sp;



%opt_params struct
opt_params.Niter = num_iter;
opt_params.gamma = gamma;

Y_avg = squeeze(mean(Y,3)); 
gpu_flag=0;

truth.A=A;
truth.source_inds=source_inds;
truth.sd=sd;
truth.sigma_m=sigma_m;
truth.q=q;

[state, params, LL_complete]= saem(init, model, opt_params, Y_avg, Y, Ns, Np, P, dist_thr, gpu_flag, whiten_flag);

res_path='/m/nbe/scratch/braintrack/elife_res';
path = res_path; 
type=3; 
fname = sprintf(strcat(path,'/state_snr_%.1f_type_%d_sim_%d.mat'),bio_snr_arr(sim),type,sim);                                                                                                                                      
save(fname, 'state');
 
fname = sprintf(strcat(path,'/params_snr_%.1f_type_%d_sim_%d.mat'),bio_snr_arr(sim),type,sim);                                                                                                                                      
save(fname, 'params');

fname = sprintf(strcat(path,'/truth_snr_%.1f_type_%d_sim_%d.mat'),bio_snr_arr(sim),type,sim);                                                                                                                                      
save(fname, 'truth');

fname = sprintf(strcat(path,'/LL_%.1f_type_%d_sim_%d.mat'),bio_snr_arr(sim),type,sim);                                                                                                                                      
save(fname, 'LL_complete');




