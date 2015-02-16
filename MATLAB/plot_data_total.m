%% Import 

cd ../../APC_Build/

data = importdata('log_iterations.txt') ;

%% Plot

L = size(data)*[1 ; 0] ;

plot(1:L , data, 'LineWidth',2)
