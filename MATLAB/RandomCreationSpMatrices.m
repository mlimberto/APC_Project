% Script to create urm and icm matrice randomly and save them in the file
% 'Matrices_examples.mat'
clear all
N_user = 50;
N_Items = 20;
N_labels = 30;
density_icm=0.01;
density_urm=0.1;

icm=sprand(N_labels,N_Items,density_icm);
[i_icm,j_icm,val_icm] = find(icm);
icm=sparse(i_icm,j_icm,ones(size(i_icm)));

urm=sprand(N_labels,N_Items,density_urm);
[i_urm,j_urm,val_urm] = find(urm);
urm=sparse(i_urm,j_urm,randi([1,5],size(i_urm)));

save('Matrices_examples.mat','icm','urm');

clear all