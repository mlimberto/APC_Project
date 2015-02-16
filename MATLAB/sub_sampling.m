%% IMPORT DATAS

clear all

[FileName,PathName] = uigetfile('*.mat','Select the MAT-file extension'); 
load(strcat(PathName,FileName)); % carica tutti i dati dal file scelto e strcat mette in fila le stringhe

%% Sottocampionamento delle labels
% Genera un vettore con tutti i tipi diversi
% di labels che abbiamo in dictionary.stemTypes (attributo della struct dictionary)

labels = unique(dictionary.stemTypes);

% decidendo noi a priori quali labels tenere, andremo a costruire un vettore
% che contiene i labels da cancellare; ad esempio teniamo solo CountryOfOrigin, Year

labels_sampling = labels([5,15]);

% sotto campionamento della matrice ICM con i labels scelti

indexes=[];

for i=1:length(labels_sampling)
    indexes=[indexes ; find(strcmp(labels_sampling(i),dictionary.stemTypes))];
end

% cancelliamo le righe della matrice le cui labels non vogliamo tenere
icm_sampling=icm;
n_rows = size(icm)*[1;0];

icm_sampling(setdiff(1:n_rows,indexes),:)=[];

new_labels = dictionary.stems(indexes);


%% Sottocampionamento del numero di utenti

N_USERS = 4000; % We want to keep only a small number of users



%% Sottocampionamento degli items
[n,m]=size(urm);
N_ITEMS =104; %% round(N_USERS*size(urm,2)/size(urm,1));

urm_sampling = urm(1:N_USERS,1:N_ITEMS);

icm_sampling = icm_sampling(:,1:N_ITEMS);

disp 'Size of URM'
size(urm_sampling)
disp 'Size of ICM'
size(icm_sampling)


%% Esportiamo ICM
% Andiamo ora a scrivere il file di testo con il sottocampionamento della
% matrice icm, in formato 'sparse'

[row_icm_sampling,column_icm_sampling,values_icm_sampling] = find(icm_sampling);
fileID_icm_sampling=fopen('icm_sampling.txt','wt');
fprintf(fileID_icm_sampling,'%d %d %d\n',[row_icm_sampling-1,column_icm_sampling-1,values_icm_sampling]');
fclose(fileID_icm_sampling);

%% Esportiamo URM

[ru,cu,vu] = find(urm_sampling);
fileID_urm_sampling=fopen('urm_sampling.txt','wt');
fprintf(fileID_urm_sampling,'%d %d %f\n',[ru-1,cu-1,vu]');
fclose(fileID_urm_sampling);
