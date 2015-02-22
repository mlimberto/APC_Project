%% Import topN recommendations

cd ../dataset/Dataset_bin_for_Tuning/

top_N = importdata('top_N.txt')

[U N] = size(top_N);

%% Import k relevant items

K = 15 ;

k_relevant = importdata('urm_bin_test.txt');

L = size(k_relevant)*[1 0]';


%% Loop on users 

precision = [];
recall = [];

for n = 1:N

local_precision = [];
local_recall = [];

for i=0:(U-1) 
    recommended = top_N(i+1 , 1:n);
    relevant = [];
    for j=1:L 
        if (k_relevant(j,1) == i )
            relevant = [relevant , k_relevant(j,2)];
        end
    end
    local_precision = [ local_precision , sum(ismember(relevant,recommended))/n ] 
    local_recall = [ local_recall , sum(ismember(relevant,recommended))/K ] 
end

precision = [precision , mean(local_precision)];
recall = [recall , mean(local_recall)];

end