function regressionanalysis

clear all; close all;
clc;
%
% Please complile the data values required in this first part of the file
% 
p=path;
curr_dir = pwd;
path(p,curr_dir);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    REQUIRED DATA TO BE FILLED BY THE USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 'fname' contains the experiment design used in the simulations
% 'dirname' is the directory where the simulated data is located
%
% dirnameHeartFailure='../results/exp_des_uniform_4_3060.HAL.280312/post_data';
% dirnameNormalVar='../results/exp_des_uniform_4_30.HAL.260312/post_data';
dirnameHeartFailure='../results/OHara2011_endo_uniform_3060/post_data';
dirnameNormalVar='../results/OHara2011_endo_uniform_30/post_data';
dirpost='../results/';
cd(dirpost);
% [sucess,message,messageid]=mkdir('regression');
[sucess,message,messageid]=mkdir('regression_Chaste');
if(sucess==0)
    fprintf('dir regression. %s. Matlab error: %s\n',message,messageid);
    return;
elseif(isempty(message))
    fprintf('dir regression. %s. Matlab error: %s\n',message,messageid);
end
%
mask_string ={'GCaL','Tauf', 'Tauf2', 'GKr', 'GKs', 'TauXs', 'GK1', 'GNaK', 'GNaCa', 'Gto', 'Jup'};
%
ranges_HF = [0.4 1; 0.4 1; 0.7 1.3; 0.7 1.3; 1 1.6; 0.4 1; 0.4 1];
ranges_NV = [0.7 1.3;0.7 1.3;0.7 1.3;0.7 1.3;0.7 1.3;0.7 1.3;0.7 1.3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get baseline values. Postprocessed values are measured as changes
% relative to the baseline.
%
cd(curr_dir);

% load HF alternans data:
cd(dirnameHeartFailure)
HF_APAlt = load('APAlt.dat');
[index_HF_APAlt_i,index_HF_APAlt_j]=find(HF_APAlt>0);
HF_CaiTAlt = load('CaAlt.dat');
[index_HF_CaAlt_i,index_HF_CaAlt_j]=find(HF_CaiTAlt>0);
cd(curr_dir);
% load NV alternans data
cd(dirnameNormalVar);
NV_APAlt = load('APAlt.dat');
[index_NV_APAlt_i,index_NV_APAlt_j]=find(NV_APAlt>0);
NV_CaiTAlt = load('CaAlt.dat');
[index_NV_CaAlt_i,index_NV_CaAlt_j]=find(NV_CaiTAlt>0);
cd(curr_dir)

% Now find the difference between calcium and AP alternans

HF_APAlt_indices = [index_HF_APAlt_i index_HF_APAlt_j];
HF_CaAlt_indices = [index_HF_CaAlt_i index_HF_CaAlt_j];
NV_APAlt_indices = [index_NV_APAlt_i index_NV_APAlt_j];
NV_CaAlt_indices = [index_NV_CaAlt_i index_NV_CaAlt_j];


cd(dirnameHeartFailure);
baseline=load('baselineBio.dat');
BCL = baseline(:,1);
clear baselineBio;

% write alternans files


% Plot alternans parameters:
%dataset_size = [16385 11];
%PlotAlternansCases(HF_APAlt_indices, HF_CaAlt_indices, NV_APAlt_indices, NV_CaAlt_indices, dataset_size, BCL)

%
% load HF parameter variations:
HFExplanatoryData=load('exp_design.dat');
mask = HFExplanatoryData(1,:);
HFExplanatoryData = HFExplanatoryData(2:end,:);

% load HF response variables
% Note, each column is for the different CLs contained in BCL
% APD80
HF_APD80 = load('APD80.dat');
clear APD80;

cd(curr_dir);
cd(dirpost);
% cd('regression');
cd('regression_Chaste');


% first write the alternans files
for i=1:length(BCL)
     write_alternans_files(i, index_HF_APAlt_i,index_HF_APAlt_j,BCL,'HFAP');
     write_alternans_files(i, index_HF_CaAlt_i,index_HF_CaAlt_j,BCL,'HFCa');
     write_alternans_files(i, index_NV_APAlt_i,index_NV_APAlt_j,BCL,'NVAP');
     write_alternans_files(i, index_NV_CaAlt_i,index_NV_CaAlt_j,BCL,'NVCa');
end


% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
[explanatory_data, biomarker] = remove_alternans(HFExplanatoryData, HF_APD80(:,i), index_HF_APAlt_i, index_HF_APAlt_j,i);
% now alternans is removed, we remove any points which had a '-1' reported
% against them in the code.
I = find(biomarker==-100);
biomarker(I) = [];
explanatory_data(I,:)=[];
% Mean centre data
explanatory_data = Normalise_Data(explanatory_data);
biomarker = Normalise_Data(biomarker);
B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
[~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
B2(:,i) = coeffs;
biomarker_cell{i} = biomarker;
explanatory_data_cell{i} = explanatory_data;

% calculate global R squared at each BCL
[num_cells,num_params] = size(explanatory_data);
coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
coeffxexplanatory = coeff_matrix.*explanatory_data;
fit = sum(coeffxexplanatory, 2);
residual = fit-biomarker;
sumsquaresresid = sum(residual.^2);
sumsquarestotal = (length(biomarker)-1)*var(biomarker);
rsquared = 1 - sumsquaresresid/sumsquarestotal;
r_squared_global(i) = rsquared;

% calculate R squared for each individual parameter.
for param = 1:length(mask)
    fit_param = coeffxexplanatory(:,param);
    residual_param = fit_param-biomarker;
    ss_residual_param = sum(residual_param.^2);
    ss_total_param = (length(biomarker)-1)*var(biomarker);
    rsquared_param = 1 - ss_residual_param/ss_total_param;
    r_squared_param(i,param) = rsquared_param;
end
    
    

end
% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('APD80regression_HF.dat','w');
fprintf(fp,'%% Regression coefficients for APD80.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

%GenerateRegressionLegend(mask, mask_string)

% %PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'HF_APD80', 'HF APD_{80}');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'HF_APD80', 'HF APD_{80}');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'HF_APD80', 'HF APD_{80}')

cd(curr_dir);
cd(dirnameHeartFailure);

% APD3080
HF_APD3080 = load('APD3080.dat');
clear APD3080;

cd(curr_dir);
cd(dirpost);
% cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(HFExplanatoryData, HF_APD3080(:,i), index_HF_APAlt_i, index_HF_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end

% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('APD3080regression_HF.dat','w');
fprintf(fp,'%% Regression coefficients for APD3080.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);


% % PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'HF_APD3080', 'HF APD_{30}/APD_{80}');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'HF_APD3080', 'HF APD_{30}/APD_{80}');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'HF_APD3080', 'HF APD_{30}/APD_{80}')
cd(curr_dir);
cd(dirnameHeartFailure);

% Cai3080
HF_Cai3080 = load('Cai3080.dat');
clear Cai3080;

cd(curr_dir);
cd(dirpost);
% cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(HFExplanatoryData, HF_Cai3080(:,i), index_HF_APAlt_i, index_HF_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
end

% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('Cai3080regression_HF.dat','w');
fprintf(fp,'%% Regression coefficients for Cai3080.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);


% % PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'HF_CaiTD3080', 'HF CaiTD_{30}/CaiTD_{80}');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'HF_CaiTD3080', 'HF CaiTD_{30}/CaiTD_{80}');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'HF_CaiTD3080', 'HF CaiTD_{30}/CaiTD_{80}')
cd(curr_dir);
cd(dirnameHeartFailure);

% CaiTD80
HF_CaiTD80 = load('CaiTD80.dat');
clear CaiTD80;

cd(curr_dir);
cd(dirpost);
% cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(HFExplanatoryData, HF_CaiTD80(:,i), index_HF_APAlt_i, index_HF_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end

% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('CaiTD80regression_HF.dat','w');
fprintf(fp,'%% Regression coefficients for CaiTD80.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

% % PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'HF_CaiTD80', 'HF CaiTD_{80}');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'HF_CaiTD80', 'HF CaiTD_{80}');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'HF_CaiTD80', 'HF CaiTD_{80}')

% Plot the distribution of the alternans parameters:
for i=1:length(BCL)
[AP_explanatory_data, ~] = get_alternans(HFExplanatoryData, HF_APD80(:,i), index_HF_APAlt_i, index_HF_APAlt_j,i);
[Ca_explanatory_data, ~] = get_alternans(HFExplanatoryData, HF_CaiTD80(:,i), index_HF_CaAlt_i, index_HF_CaAlt_j,i);
if (isempty(AP_explanatory_data)==0 || isempty(Ca_explanatory_data)==0)
    % PlotAlternansDistributions(AP_explanatory_data, Ca_explanatory_data, ranges_HF, mask, mask_string,BCL(i), 'HF', 'HF Alternans')
end   
end

cd(curr_dir);
cd(dirnameHeartFailure);

% DAPCaiT
HF_DAPCaiT = load('DAPCaiT.dat');
clear DAPCaiT;

cd(curr_dir);
cd(dirpost);
% cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(HFExplanatoryData, HF_DAPCaiT(:,i), index_HF_APAlt_i, index_HF_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end

% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('DAPCaiTregression_HF.dat','w');
fprintf(fp,'%% Regression coefficients for DAPCaiT.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

%PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'HF_DAPCaiT', 'HF DAPCaiT');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'HF_DAPCaiT', 'HF DAPCaiT');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'HF_DAPCaiT', 'HF DAPCaiT')


cd(curr_dir);
cd(dirnameHeartFailure);

% SysCai
HF_SysCai = load('SysCai.dat');
clear SysCai;

cd(curr_dir);
cd(dirpost);
%cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(HFExplanatoryData, HF_SysCai(:,i), index_HF_APAlt_i, index_HF_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end

% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('SysCairegression_HF.dat','w');
fprintf(fp,'%% Regression coefficients for SysCai.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

%PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'HF_CaiTmax', 'HF CaiT_{max}');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'HF_CaiTmax', 'HF CaiT_{max}');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'HF_CaiTmax', 'HF CaiT_{max}')

cd(curr_dir);

% 
% NOW REPEAT FOR NORMAL VARIATION.
%

cd(dirnameNormalVar);
baseline=load('baselineBio.dat');
BCL = baseline(:,1);
clear baselineBio;
%
% load NV parameter variations:
NVExplanatoryData=load('exp_design.dat');
mask = NVExplanatoryData(1,:);
NVExplanatoryData = NVExplanatoryData(2:end,:);

% load NV response variables
% Note, each column is for the different CLs contained in BCL
% APD80
NV_APD80 = load('APD80.dat');
clear APD80;

cd(curr_dir);
cd(dirpost);
%cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(NVExplanatoryData, NV_APD80(:,i), index_NV_APAlt_i, index_NV_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end
% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('APD80regression_NV.dat','w');
fprintf(fp,'%% Regression coefficients for APD80.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

%PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'NV_APD80', 'NV APD_{80}');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'NV_APD80', 'NV APD_{80}');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'NV_APD80', 'NV APD_{80}')

cd(curr_dir);
cd(dirnameNormalVar);

% APD3080
NV_APD3080 = load('APD3080.dat');
clear APD3080;

cd(curr_dir);
cd(dirpost);
%cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(NVExplanatoryData, NV_APD3080(:,i), index_NV_APAlt_i, index_NV_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end
% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('APD3080regression_NV.dat','w');
fprintf(fp,'%% Regression coefficients for APD3080.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

%PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'NV_APD3080', 'NV APD_{30}/APD_{80}');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'NV_APD3080', 'NV APD_{30}/APD_{80}');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'NV_APD3080', 'NV APD_{30}/APD_{80}')
%

cd(curr_dir);
cd(dirnameNormalVar);

% Cai3080
NV_Cai3080 = load('Cai3080.dat');
clear Cai3080;

cd(curr_dir);
cd(dirpost);
% cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(NVExplanatoryData, NV_Cai3080(:,i), index_NV_APAlt_i, index_NV_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end

% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('Cai3080regression_NV.dat','w');
fprintf(fp,'%% Regression coefficients for Cai3080.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

%PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'NV_CaiTD3080', 'NV CaiTD3080');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'NV_CaiTD3080', 'NV CaiTD3080');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'NV_CaiTD3080', 'NV CaiTD3080')


cd(curr_dir);
cd(dirnameNormalVar);

% CaiTD80
NV_CaiTD80 = load('CaiTD80.dat');
clear CaiTD80;

cd(curr_dir);
cd(dirpost);
% cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(NVExplanatoryData, NV_CaiTD80(:,i), index_NV_APAlt_i, index_NV_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end
% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('CaiTD80regression_NV.dat','w');
fprintf(fp,'%% Regression coefficients for CaiTD80.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

%PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'NV_CaiTD80', 'NV CaiTD_{80}');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'NV_CaiTD80', 'NV CaiTD_{80}');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'NV_CaiTD80', 'NV CaiTD_{80}')

% For CaiTD80, we now plot the distribution of the alternans parameters:
for i=1:length(BCL)
[AP_explanatory_data, ~] = get_alternans(NVExplanatoryData, NV_APD80(:,i), index_NV_APAlt_i, index_NV_APAlt_j,i);
[Ca_explanatory_data, ~] = get_alternans(NVExplanatoryData, NV_CaiTD80(:,i), index_NV_CaAlt_i, index_NV_CaAlt_j,i);
if (isempty(AP_explanatory_data)==0 || isempty(Ca_explanatory_data)==0)
   % PlotAlternansDistributions(AP_explanatory_data, Ca_explanatory_data, ranges_NV, mask, mask_string,BCL(i), 'NV', 'NV Alternans')
end   
end

cd(curr_dir);
cd(dirnameNormalVar);

% DAPCaiT
NV_DAPCaiT = load('DAPCaiT.dat');
clear DAPCaiT;

cd(curr_dir);
cd(dirpost);
% cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(NVExplanatoryData, NV_DAPCaiT(:,i), index_NV_APAlt_i, index_NV_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end

% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('DAPCaiTregression_NV.dat','w');
fprintf(fp,'%% Regression coefficients for DAPCaiT.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

%PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'NV_DAPCaiT', 'NV DAPCaiT');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'NV_DAPCaiT', 'NV DAPCaiT');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'NV_DAPCaiT', 'NV DAPCaiT')

cd(curr_dir)
cd(dirnameNormalVar);

% SysCai
NV_SysCai = load('SysCai.dat');
clear SysCai;

cd(curr_dir);
cd(dirpost);
%cd('regression');
cd('regression_Chaste');

% calculate in two different ways. Should be equal.
B1 = zeros(length(mask), length(BCL));
B2 = zeros(length(mask)+1, length(BCL));

r_squared_global = zeros(length(BCL),1);
r_squared_param = zeros(length(BCL),length(mask));

for i=1:length(BCL)
    [explanatory_data, biomarker] = remove_alternans(NVExplanatoryData, NV_SysCai(:,i), index_NV_APAlt_i, index_NV_APAlt_j,i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I = find(biomarker==-100);
    biomarker(I) = [];
    explanatory_data(I,:)=[];
    % Mean centre data
    explanatory_data = Normalise_Data(explanatory_data);
    biomarker = Normalise_Data(biomarker);
    B1(:,i) = (explanatory_data'*explanatory_data)\explanatory_data'*biomarker;
    [~,~,~,~,coeffs] = plsregress(explanatory_data,biomarker, 7);
    B2(:,i) = coeffs;
    biomarker_cell{i} = biomarker;
    explanatory_data_cell{i} = explanatory_data;
    
    % calculate global R squared at each BCL
    [num_cells,num_params] = size(explanatory_data);
    coeff_matrix = ones(num_cells, 1)*coeffs(2:end)';
    coeffxexplanatory = coeff_matrix.*explanatory_data;
    fit = sum(coeffxexplanatory, 2);
    residual = fit-biomarker;
    sumsquaresresid = sum(residual.^2);
    sumsquarestotal = (length(biomarker)-1)*var(biomarker);
    rsquared = 1 - sumsquaresresid/sumsquarestotal;
    r_squared_global(i) = rsquared;
    
    % calculate R squared for each individual parameter.
    for param = 1:length(mask)
        fit_param = coeffxexplanatory(:,param);
        residual_param = fit_param-biomarker;
        ss_residual_param = sum(residual_param.^2);
        ss_total_param = (length(biomarker)-1)*var(biomarker);
        rsquared_param = 1 - ss_residual_param/ss_total_param;
        r_squared_param(i,param) = rsquared_param;
    end
    
end

% Check these answers are equal:
assert(max(max(B1-B2(2:end,:))) < 1e-12);

% transpose B1 so that each row is the coefficients for a cycle length:
coefficients = B1';

% write results into a data file
fp=fopen('SysCairegression_NV.dat','w');
fprintf(fp,'%% Regression coefficients for SysCai.\n');
fprintf(fp,'%% Explanatory and response data are both mean centred and normalised by standard deviation\n');
fprintf(fp,'%% 1st row BCL\n');
fprintf(fp,'%% 2nd row mask\n');
fprintf(fp,'%% Subsequent rows: Coefficients corresponding to mask at each BCL.\n');
fprintf(fp,'%6.0f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%6.0f ',mask);
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',coefficients(i,:));
    fprintf(fp,'\n');
end

for j=1:length(mask)
    fprintf(fp,'%1.3f & ',coefficients(:,j)');
    fprintf(fp,'\\ \n');
end

fprintf(fp, '%% global R squared coefficients at each BCL\n');
fprintf(fp, '%12.4e ', r_squared_global');
fprintf(fp,'\n');
fprintf(fp, '%% individual R squared coefficients at each BCL\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',r_squared_param(i,:));
    fprintf(fp,'\n');
end
fclose(fp);

%PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, 'NV_CaiTmax', 'NV CaiT_{max}');
PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, 'NV_CaiTmax', 'NV CaiT_{max}');
%PlotOutputVsBiomarker(BCL, coefficients, biomarker_cell, explanatory_data_cell, mask, mask_string, 'NV_CaiTmax', 'NV CaiT_{max}')

cd(curr_dir)

end

function PlotCoefficientBarCharts(BCL, coefficients, mask, mask_string, FileNameString, FigTitleString)
% Plot bar charts for each cycle length
for i=1:length(BCL)
    fig=figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    % make bar chart
    bar(coefficients(i,:));
    % Set the correct labels for the bar chart
    set(gca,'XTickLabel',mask_string(mask+1));
    % Add title
    title_text=sprintf('%s, BCL=%.0f ms',FigTitleString,BCL(i));
    title(title_text);
    % remove box
    set(gca,'Box','off');
    % set font size for all text
    hAll = findall(gcf);
    for idx = 1 : length(hAll)
        try
            set(hAll(idx),'fontsize',30);
        catch
            % nothing
        end
    end
    % Set width of axes
    set(gca, 'linewidth', 2);
    % make background transparent
    set(gca,'color','none');
    % save to file using export_fig
    sname=sprintf('regression_%s_BCL%.0f',FileNameString,BCL(i));
    export_fig([sname '.png'], '-png', '-transparent');
    export_fig([sname '.eps'], '-eps', '-transparent');
    close all
end

end

function GenerateRegressionLegend(mask, mask_string)
for i=1:4
    fig = figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    if i==1
        legend_index = [1 2 5];
    elseif i==2
        legend_index = [1 2 4 5 7];
    elseif i==3
        legend_index = [1 2 7];
    elseif i==4
        legend_index = [1 7];
    else
        legend_index =[];
    end
    for j=1:7;
        BCL_plot{j} = [];
        coefficients_plot{j} = [];
    end
    for j=1:length(legend_index)
        BCL_plot{legend_index(j)} = 1:10;
        coefficients_plot{legend_index(j)} = ones(1,10);
    end
    hold on
    plot(BCL_plot{1}, coefficients_plot{1}, 'kx', BCL_plot{2}, coefficients_plot{2}, 'ko',...
        BCL_plot{3}, coefficients_plot{3}, 'k+',BCL_plot{4}, coefficients_plot{4}, 'k*',...
        BCL_plot{5}, coefficients_plot{5}, 'ks',BCL_plot{6}, coefficients_plot{6}, 'kd',...
        BCL_plot{7}, coefficients_plot{7}, 'k^', 'MarkerSize', 30, 'LineWidth', 4);
    hold off
    set(gca, 'box', 'off');
    
    set(gca, 'ylim', [0 10]);
    
    regression_leg=legend(mask_string{mask(legend_index)+1});
    set(regression_leg, 'box', 'off');
    hAll = findall(gcf);
    for idx = 1 : length(hAll)
        try
            set(hAll(idx),'fontsize',40);
        catch
            % nothing
        end
    end
    sname=sprintf('regression_legend_%0.f', i);
    export_fig([sname '.png'], '-png', '-transparent');
    export_fig([sname '.eps'], '-eps', '-transparent');
    close all
end

end

function PlotCoefficientsAllBCLs(BCL, coefficients, mask, mask_string, FileNameString, FigTitleString)
% Plot graphs for all cycle lengths:

fig=figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
% We just plot the coefficients with values > 0.2.
legend_index = [];
for i=1:length(mask)
    index = find(abs(coefficients(:,i)) > 0.2);
    if ~isempty(index)
        BCL_plot{i} = BCL;
        coefficients_plot{i} = coefficients(:,i);
    else
        BCL_plot{i} = [];
        coefficients_plot{i} = [];
    end
    %BCL_plot{i} = BCL(index);
    %coefficients_plot{i} = coefficients(index,i);
    if ~isempty(coefficients_plot{i})
        legend_index = [legend_index i];
    end
end
% add rectangle

rect=rectangle('Position',[305,-0.2,1195 0.4],'FaceColor', 0.7*[1 1 1], 'EdgeColor', 0.7*[1 1 1]);
hold on
plot(BCL_plot{1}, coefficients_plot{1}, 'kx', BCL_plot{2}, coefficients_plot{2}, 'ko',...
    BCL_plot{3}, coefficients_plot{3}, 'k+',BCL_plot{4}, coefficients_plot{4}, 'k*',...
    BCL_plot{5}, coefficients_plot{5}, 'ks',BCL_plot{6}, coefficients_plot{6}, 'kd',...
    BCL_plot{7}, coefficients_plot{7}, 'k^', 'MarkerSize', 50, 'LineWidth', 8);
hold off
%legend_string = mask_string{mask(legend_index)+1};
%regression_leg=legend(mask_string{mask(legend_index)+1});
%set(regression_leg, 'Position', get(regression_leg, 'Position') + [0.05 0 0 0]);
hold on
plot(BCL_plot{1}, coefficients_plot{1}, 'k-', BCL_plot{2}, coefficients_plot{2}, 'k-',...
    BCL_plot{3}, coefficients_plot{3}, 'k-',BCL_plot{4}, coefficients_plot{4}, 'k-',...
    BCL_plot{5}, coefficients_plot{5}, 'k-',BCL_plot{6}, coefficients_plot{6}, 'k-',...
    BCL_plot{7}, coefficients_plot{7}, 'k-', 'LineWidth', 8);
hold off
% Add title
%title_text=sprintf('Regression Coefficients, %s', FigTitleString);
%title(title_text);
% set axis length
ylimits = get(gca,'ylim');
axis([300 1500 -1 1]);
% set axis labels
xlabel('BCL (ms)')
% remove box
set(gca,'Box','off');

% set font size for all text
hAll = findall(gcf);
for idx = 1 : length(hAll)
    try
        set(hAll(idx),'fontsize',60);
    catch
        % nothing
    end
end
% Set width of axes
set(gca, 'linewidth', 6);
% make background transparent
set(gca,'color','none');
% save to file using export_fig
sname=sprintf('regression_allBCLs_%s', FileNameString);
export_fig([sname '.png'], '-png', '-transparent');
export_fig([sname '.eps'], '-eps', '-transparent');
close all
end

function PlotOutputVsBiomarker(BCL, coefficients, biomarker_normalised, input_normalised, mask, mask_string, FileNameString, FigTitleString)
for i=1:length(BCL)
    biomarker = biomarker_normalised{i};
    input_BCL = input_normalised{i};
    for j=1:length(mask)
        fig=figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
        input = squeeze(input_BCL(:,j));
        plot(input, biomarker, 'k.', 'MarkerSize', 6, 'LineWidth', 4);
        limits=get(gca, 'xlim');
        xdata =linspace(limits(1),limits(2),15);
        regressdata = coefficients(i,j)*xdata;
        hold on
        plot(xdata, regressdata,'color', 0.3*[1 1 1], 'LineStyle', '-', 'LineWidth', 6);
        plot(xdata, regressdata,'color', 0.6*[1 1 1], 'LineStyle', '--','LineWidth', 6);
        hold off
        %legend('data','Regression');
        % Add title
        %title_text=sprintf('%s vs %s, BCL=%.0f ms',FigTitleString,mask_string{mask(j)+1},BCL(i));
        %title(title_text);
        % remove box
        set(gca,'Box','off');
        % set font size for all text
        
        xlabel(sprintf( '%s (normalised)', mask_string{mask(j)+1}));
        ylabel(sprintf('%s (normalised)', FigTitleString));
        
        hAll = findall(gcf);
        for idx = 1 : length(hAll)
            try
                set(hAll(idx),'fontsize',40);
            catch
                % nothing
            end
        end
        % Set width of axes
        set(gca, 'linewidth',3);
        % make background transparent
        set(gca,'color','none');
        % save to file using export_fig
        sname=sprintf('regression_%s%s_BCL%.0f',FileNameString,mask_string{mask(j)+1},BCL(i));
        export_fig([sname '.png'], '-png', '-transparent');
        export_fig([sname '.eps'], '-eps', '-transparent');
        close all
    end
end

end

function write_alternans_files(current_CL_index, index_alt_i,index_alt_j,BCLs,FileNameString)
if any(index_alt_j == current_CL_index)
    I=find(index_alt_j==current_CL_index);
    alternans_cases=index_alt_i(I);
    filename = [FileNameString '_AlternansCases_' num2str(BCLs(current_CL_index)) '.dat'];
    fp=fopen(filename,'w');
    fprintf(fp,'%% Experiment_Numbers\n');
    fprintf(fp, '%5.0f \n', alternans_cases);
    fclose(fp);
end
end

function [explanatory_noalternans, biomarker_noalternans] = remove_alternans(explanatory, biomarker, index_alt_i, index_alt_j, current_CL_index)
explanatory_noalternans=explanatory;
biomarker_noalternans=biomarker;
if any(index_alt_j == current_CL_index)
    I=find(index_alt_j==current_CL_index);
    explanatory_noalternans(index_alt_i(I),:) = [];
    biomarker_noalternans(index_alt_i(I))=[];
end
end

function [explanatory_alternans, biomarker_alternans] = get_alternans(explanatory, biomarker, index_alt_i, index_alt_j, current_CL_index)
explanatory_alternans=[];
biomarker_alternans=[];

if any(index_alt_j == current_CL_index)
    I=find(index_alt_j==current_CL_index);
    explanatory_alternans = explanatory(index_alt_i(I),:);
    biomarker_alternans= biomarker(index_alt_i(I));
end
end

function PlotAlternansDistributions(AP_explanatory, Ca_explanatory,ranges, mask, mask_string,cycle_length, FileNameString, FigTitleString)
figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);


size(AP_explanatory)
size(Ca_explanatory)

APandCaiTAlternans = intersect(AP_explanatory, Ca_explanatory, 'rows');
APnotCaiTAlternans = setdiff(AP_explanatory, Ca_explanatory, 'rows');
CaiTnotAPAlternans = setdiff(Ca_explanatory, AP_explanatory, 'rows');

AllAlternans = union(Ca_explanatory, AP_explanatory,'rows');

[num_cases,~] = size(AllAlternans);

% ensure number of intervals is divisible by 4
num_bins = max(round(sqrt(num_cases)), 10);

num_bins = num_bins - mod(num_bins, 4);

num_bins = 8;

% write results into a data file
filename = [FileNameString '_AlternansHistograms_' num2str(cycle_length) '.dat'];
fp=fopen(filename,'w');
fprintf(fp,'%% Histogram values\n');
fprintf(fp,'%% Number of bins\n');
fprintf(fp, '%6.0f\n', num_bins);
fprintf(fp,'%% For each parameter:\n');
fprintf(fp,'%% 1st row AP and CaiT alternans\n');
fprintf(fp,'%% 2nd row AP alternans only\n');
fprintf(fp,'%% 3rd row Ca alternans only\n');

% get max value of histogram:

max_vector = zeros(1, num_bins+1);
for i=1:length(mask)
    edges = linspace(min(ranges(i,1), ranges(i,2)),max(ranges(i,1), ranges(i,2)),num_bins+1);
    Hist_i = histc(AllAlternans(:,i),edges);
    max_vector(i) = max(Hist_i);
end
hist_y_limit = max(max_vector);

for i=1:length(mask)
    % Get correct bin spacing:
   
    edges = linspace(min(ranges(i,1), ranges(i,2)),max(ranges(i,1), ranges(i,2)),num_bins+1);
    
    % this is an integer as num_bins is divisible by 4;
    put_label_here = num_bins/4;
    
    for j=0:length(edges)-1
        if mod(j,put_label_here) == 0
            ticks{j+1} = num2str(edges(j+1));
        else
            ticks{j+1} = '';
        end
    end
    
    [APandCaiT,~] = histc(APandCaiTAlternans(:,i),edges);
    [APnotCaiT,~] = histc(APnotCaiTAlternans(:,i),edges);
    [CaiTnotAP,~] = histc(CaiTnotAPAlternans(:,i),edges);
       
    if(size(APandCaiT)~= [num_bins 1])
        APandCaiT = APandCaiT';
    end
    if(size(APnotCaiT)~= [num_bins 1])
        APnotCaiT = APnotCaiT';
    end
    if(size(CaiTnotAP)~= [num_bins 1])
        CaiTnotAP = CaiTnotAP';
    end
    
    bardata = [APandCaiT, APnotCaiT, CaiTnotAP];

    subplot(4,2,i);bar(bardata,1,'stacked');
        
    set(gca, 'ylim', [0 hist_y_limit]);
    set(gca,'xtickmode', 'auto');
    set(gca,'xtick',(0.5:length(edges)-0.5));
    set(gca, 'xticklabel', ticks);

    xlabel([mask_string{mask(i)+1} ' scale']);
    % remove box
    set(gca,'Box','off');
    % Set width of axes
    set(gca, 'linewidth', 2);
    % make background transparent
    set(gca,'color','none');
    % correct axes lengths
    set(gca, 'xlim', [0.5 length(edges)-0.5]);
    set(gca, 'TickDir', 'out');
    
    greyscale_colourmap = [ 0.2*[1 1 1]; 0.5*[1 1 1]; 0.7*[1 1 1]]; % both, AP, Ca'
    colormap(greyscale_colourmap);
   
    clear ticks
    
    fprintf(fp, mask_string{mask(i)+1});
    fprintf(fp, '\n');
    fprintf(fp, '%6.0f &', APandCaiT);
    fprintf(fp, '\n');
    fprintf(fp, '%6.0f &', APnotCaiT);
    fprintf(fp, '\n');
    fprintf(fp, '%6.0f &', CaiTnotAP);
    fprintf(fp, '\n');
    total_alternans = sum(APandCaiT)+sum(APnotCaiT)+sum(CaiTnotAP);
    fprintf(fp, 'num cases %6.0f', total_alternans);
    fprintf(fp, '\n');
    fprintf(fp, '%1.5f ', APandCaiT/total_alternans);
    fprintf(fp, '\n');
    fprintf(fp, '%1.5f ', APnotCaiT/total_alternans);
    fprintf(fp, '\n');
    fprintf(fp, '%1.5f ', CaiTnotAP/total_alternans);
    fprintf(fp, '\n');
    
    
end

fclose(fp);

%figtitle = sprintf('%s, BCL %.0f ms',FigTitleString, cycle_length);
%[~,sup_handle]=suplabel(figtitle, 't');

% set font size for all text
hAll = findall(gcf);
for idx = 1 : length(hAll)
    try
        set(hAll(idx),'fontsize',30);
    catch
        % nothing
    end
end

%set(sup_handle, 'Position', get(sup_handle, 'Position')+[0 -0.05 0]);

% save to file using export_fig
sname=sprintf('Alternans_%s_BCL=%.0f', FileNameString, cycle_length);
%export_fig([sname '.png'], '-png', '-transparent');
%export_fig([sname '.eps'], '-eps', '-transparent');
close all

end

function normalised_data = Normalise_Data(data)
[data_rows, ~] = size(data);
normalised_data = (data - (ones(data_rows,1)*mean(data)))./(ones(data_rows,1)*sqrt(var(data)));
end

% function PlotAlternansCases(HF_APAlt_indices, HF_CaAlt_indices, NV_APAlt_indices, NV_CaAlt_indices, dataset_size, BCL)
% 
% fig=figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
% 
% 
% 
% % 
% % HF_APnotCaiTAlternans = setdiff(HF_APAlt_indices, HF_CaAlt_indices, 'rows');
% % HF_CaiTnotAPAlternans = setdiff(HF_CaAlt_indices, HF_APAlt_indices, 'rows');
% % NV_APnotCaiTAlternans = setdiff(NV_APAlt_indices, NV_CaAlt_indices, 'rows');
% % NV_CaiTnotAPAlternans = setdiff(NV_CaAlt_indices, NV_APAlt_indices, 'rows');
% % 
% % HF_APandCaiTAlternans = intersect(HF_APAlt_indices, HF_CaAlt_indices, 'rows');
% % NV_APandCaiTAlternans = intersect(NV_APAlt_indices, NV_CaAlt_indices, 'rows');
% % 
% % 
% % fig=figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
% % subplot(1,2,1); 
% % plot(NV_APandCaiTAlternans(:,2),NV_APandCaiTAlternans(:,1), 'ko',NV_APnotCaiTAlternans(:,2),...
% %     NV_APnotCaiTAlternans(:,1), 'rx', NV_CaiTnotAPAlternans(:,2), NV_CaiTnotAPAlternans(:,1), 'b+', 'LineWidth', 2, 'MarkerSize',1);
% % legend('Both', 'AP Alternans', 'CaiT Alternans');
% % title('Alternans, Non-Failing cases')
% % subplot(1,2,2); 
% % plot(HF_APandCaiTAlternans(:,2),HF_APandCaiTAlternans(:,1), 'ko',HF_APnotCaiTAlternans(:,2),...
% %     HF_APnotCaiTAlternans(:,1), 'rx', HF_CaiTnotAPAlternans(:,2), HF_CaiTnotAPAlternans(:,1), 'b+', 'LineWidth', 2, 'MarkerSize',1);
% % legend('Both', 'AP Alternans', 'CaiT Alternans');
% % title('Alternans, Failing cases')
% end
