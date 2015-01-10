%% Script per estrarre un test set dalla matrice URM

% Load .mat file
[FileName,PathName] = uigetfile('*.mat','Select the MAT-file where the matrices urm and icm are stored'); 
load(strcat(PathName,FileName));

%% Sample the data randomly

[i_urm,j_urm,val_urm] = find(urm);

N_INPUTS = size(i_urm)*[1;0];

indexes = datasample(1:N_INPUTS,idivide(int32(N_INPUTS),10,'floor'),'Replace',false);

i_urm_test = i_urm(indexes);
i_urm(indexes) = [];

j_urm_test = j_urm(indexes);
j_urm(indexes) = [];

val_urm_test = val_urm(indexes);
val_urm(indexes) = [];


%% Export the matrices

% URM_TRAINING
fileID_urm=fopen('urm_converted.txt','wt');
fprintf(fileID_urm,'%d %d %f\n',[i_urm-1,j_urm-1,val_urm]');
fclose(fileID_urm);

% ICM
[i_icm,j_icm,values_icm] = find(icm);
fileID_icm=fopen('icm_converted.txt','wt');
fprintf(fileID_icm,'%d %d %d\n',[i_icm-1,j_icm-1,values_icm]');
fclose(fileID_icm);

% URM_TEST
fileID_urmt=fopen('urm_test.txt','wt');
fprintf(fileID_urmt,'%d %d %f\n',[i_urm_test-1,j_urm_test-1,val_urm_test]');
fclose(fileID_urmt);

