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


%Suddivisione di urm_tuning_sampling in urm_tuning_tr e urm_tuning_val togliendo a caso il 10% degli elementi della matrice

urm_bin_tuning_tr = urm_sampling
vec = find(urm_bin_tuning_tr);
f = 0.1; % desired fraction
n = round(f*length(vec));
vec_change = randsample(vec,n);
urm_bin_tuning_tr(vec_change) = 0;
urm_bin_tuning_val=urm_sampling-urm_bin_tuning_tr;


%% Esportiamo URM_tuning_tr

[ru,cu,vu] = find(urm_bin_tuning_tr);
fileID_urm_tr=fopen('urm_bin_tuning_tr.txt','wt');
fprintf(fileID_urm_tr,'%d %d %f\n',[ru-1,cu-1,vu]');
fclose(fileID_urm_tr);

%% Esportiamo URM_tuning_val

[ru,cu,vu] = find(urm_bin_tuning_val);
fileID_urm_test=fopen('urm_bin_tuning_val.txt','wt');
fprintf(fileID_urm_test,'%d %d %f\n',[ru-1,cu-1,vu]');
fclose(fileID_urm_test);


