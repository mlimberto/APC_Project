
% IMPORT DATAS

clear all

[FileName,PathName] = uigetfile('*.mat','Select the MAT-file extension'); 
load(strcat(PathName,FileName)); % carica tutti i dati dal file scelto e strcat mette in fila le stringhe

% QUALI E QUANTI SONO LE LABELS? Genera un vettore con tutti i tipi diversi
% di labels che abbiamo in dictionary.stemTypes (attributo della struct dictionary)

labels = unique(dictionary.stemTypes);

% decidendo noi a priori quali labels tenere, andremo a costruire un vettore
% che contiene i labels da cancellare; ad esempio teniamo solo ShowType, TitleFull, Year

labels_sampling = labels( 1 : end - 3 );

% sotto campionamento della matrice ICM con i labels scelti

indexes=[];

for i=1:length(labels_sampling)
    indexes=[indexes ; find(strcmp(labels_sampling(i),dictionary.stemTypes))];
end
icm_sampling=icm;

% cancelliamo le righe della matrice le cui labels non vogliamo tenere

icm_sampling(indexes,:)=[];

% Andiamo ora a scrivere il file di testo con il sottocampionamento della
% matrice icm, in formato 'sparse'

[row_icm_sampling,column_icm_sampling,values_icm_sampling] = find(icm_sampling);
fileID_icm_sampling=fopen('icm_sampling.txt','wt');
fprintf(fileID_icm_sampling,'%d %d %d\n',[row_icm_sampling-1,column_icm_sampling-1,values_icm_sampling]');
fclose(fileID_icm_sampling);
