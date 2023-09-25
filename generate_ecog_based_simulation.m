function [Y, G_mag, sigma_meas] = run_ecog_based_simulation(sub_no, meas_snr, bio_snr)
ecog_data_path = '/opt/lintula/worktmp/nbnapu/ecog_sim_data';
addpath(strcat(pwd, '/EcoG_simu'));
addpath(strcat(pwd,'/external/mvar_functions'))
addpath(strcat(pwd, '/external/mvar_functions/external/arfit'))
addpath(strcat(pwd,'/misc'))
addpath(ecog_data_path)
addpath(strcat(pwd,'/MVAR_simu'))
addpath(strcat(pwd,'/two_step_est/source_localizers'))




%%load data
if sub_no < 10
    load(sprintf(strcat(ecog_data_path,'/sub_0%d_model.mat'), sub_no))
else
    load(sprintf(strcat(ecog_data_path,'/sub_00%d_model.mat'), sub_no))
end
    
load(sprintf(strcat(pwd, '/EcoG_simu/ERP_face_bl_resamp.mat')))


mesh=mm.mesh{1,1};

anat_no = [24 14 62]; %LOC-FG-STS
mesh.p=mesh.p*1000; %in mm


%select sources
fg_pnts = mesh.p(included_sp(anat_ds.parcels{anat_no(2),1}),:);
loc_pnts = mesh.p(included_sp(anat_ds.parcels{anat_no(1),1}),:);
sts_pnts = mesh.p(included_sp(anat_ds.parcels{anat_no(3),1}),:);

cent=mean(loc_pnts);
[~, id_loc] = find_nearest_id(cent, loc_pnts(:,:));

cent=mean(fg_pnts);
[~, id_fg] = find_nearest_id(cent, fg_pnts(:,:));

cent=mean(sts_pnts);
[~, id_sts] = find_nearest_id(cent, sts_pnts(:,:));

lf_id(1) = anat_ds.parcels{24,1}(id_loc(1));
lf_id(2) = anat_ds.parcels{14,1}(id_fg(1));
lf_id(3) = anat_ds.parcels{62,1}(id_sts(1));



mag=1:3:306;
G_mag=G(mag,:);

x_trials_all=permute(ERP_face_bl.trial, [2 3 1]);
elec = [77 73 62];
t_beg=0;
t_end=0.5;
t_diff=abs(ERP_face_bl.time-t_beg);
id_beg=find(t_diff==min(t_diff));
t_diff=abs(ERP_face_bl.time-t_end);
id_end=find(t_diff==min(t_diff));
x = x_trials_all(elec,id_beg:id_end,:)*10^-9;

N=size(G_mag,1);
J=size(x,3);
T=size(x,2);
Y=zeros(N,T,J);
%add brain and measurement noise
for j=1:1:J
    Y_nn = G_mag(:,lf_id)*x(:,:,j);
    [Ybn, ~] = add_brain_noise(Y_nn, G_mag, bio_snr, 2000);
    [Y(:,:,j), sigma_meas, ~] = add_meas_noise(Ybn, meas_snr);
end


end

