%% Matlab script to convert a Sparse Matrix from the .mat format into txt.
% 
%%
function SpMat_ConvertWithLocations

% INPUT
[FileName,PathName] = uigetfile('*.mat','Select the MAT-file where the matrices urm and icm are stored'); 
load(strcat(PathName,FileName));

% check the existence of the icm and urm matrices
if ~exist('urm','var')
    disp ('ERROR : The matrix urm cannot be found');
elseif (~exist('icm','var'))
    disp ('ERROR : The matrix icm cannot be found');
else
    % WRITE THE .txt
    [i_urm,j_urm,values_urm] = find(urm);
    fileID_urm=fopen('urm_converted.txt','wt');
    fprintf(fileID_urm,'%d %d %f\n',[i_urm-1,j_urm-1,values_urm]');
    fclose(fileID_urm);

    [i_icm,j_icm,values_icm] = find(icm);
    fileID_icm=fopen('icm_converted.txt','wt');
    fprintf(fileID_icm,'%d %d %d\n',[i_icm-1,j_icm-1,values_icm]');
    fclose(fileID_icm);
end

end