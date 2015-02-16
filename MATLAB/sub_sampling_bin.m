%% IMPORT DATAS

clear all
close all

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
icm_unique=icm;
n_rows = size(icm)*[1;0];

icm_unique(setdiff(1:n_rows,indexes),:)=[];

new_labels = dictionary.stems(indexes);


%% PARAMETRI
relevant_Val=3; %% valutazione minima rilevante
M=15; %numero di valutazioni da togliere ad ogni utente nella urm_sampling
N_USERS=4000;
N_ITEMS=1000; %%round(N_USERS*size(urm,2)/size(urm,1));

%% BINARIZZAZIONE urm
urm_bin=(urm>=relevant_Val);    %consideriamo solo le valutazioni maggiori o uguali di relevant_Val

%% Riordinamento degli utenti in ordine decrescente di numero votazioni rilevanti
[r_urm_bin,c_urm_bin,v_urm_bin] = find(urm_bin);
out_bin=[(1:size(urm_bin,1))',histc(r_urm_bin,(1:size(urm_bin,1))')]; %matrice che su ogni riga i ha [utente(i),n_votationi dell'utente(i)]
[n_val_sorted_bin,index_sorted_bin]=sort(out_bin(:,2),'descend');

%% Estrazione sottomatrici campionate
index_users=sort(index_sorted_bin(1:N_USERS));
perm_index_items=randperm(size(urm_bin,2));
index_items=perm_index_items(1:N_ITEMS);
urm_sampling = urm_bin(index_users,index_items);

icm_sampling = icm_unique(:,index_items);

disp 'Size of URM'
size(urm_sampling)
disp 'Size of ICM'
size(icm_sampling)
%Suddivisione di urm_sampling in urm_tr e urm_test

%Per prima cosa, controlliamo che per ogni utente di urm_sampling abbiamo
%almeno M valutazioni rilevanti
[r_urm_s,c_urm_s,v_urm_s] = find(urm_sampling);
out_s=[(1:size(urm_sampling,1))',histc(r_urm_s,(1:size(urm_sampling,1))')];

if min(out_s(:,2))<M
    disp('WARNING : for some users m is grater than the number of relevant evaluations');
end
%%
urm_tr=urm_sampling;
for i=1:size(urm_sampling,1)
    clear i_nnz j_nnz v_nnz
    [i_nnz,j_nnz,v_nnz]=find(urm_sampling(i,:));
    perm_index_nnz=randperm(length(j_nnz));
    
    urm_tr(i,j_nnz(perm_index_nnz(1:20)))=0;
end

urm_test=urm_sampling-urm_tr;

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

%% Esportiamo URM_tr

[ru,cu,vu] = find(urm_tr);
fileID_urm_tr=fopen('urm_tr.txt','wt');
fprintf(fileID_urm_tr,'%d %d %f\n',[ru-1,cu-1,vu]');
fclose(fileID_urm_tr);

%% Esportiamo URM_test

[ru,cu,vu] = find(urm_test);
fileID_urm_test=fopen('urm_test.txt','wt');
fprintf(fileID_urm_test,'%d %d %f\n',[ru-1,cu-1,vu]');
fclose(fileID_urm_test);
