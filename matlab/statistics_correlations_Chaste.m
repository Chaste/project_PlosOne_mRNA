%
% Computes the statistics and histograms for the following biomarkers:
% APD80; AP triangulation as ADP80-ADP30; CaiTD;
% Cai triangulation as CaiTD80-CaiTD30; AP Cai delay
%
% Data is normalized with respect to the baseline value
%
% postprocessing program
%
function statistics_correlations
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
export_figures = true;

% dirnameNatVar='../results/exp_des_uniform_4_30.HAL.260312/post_data';
% dirnameHFNoVar='../results/exp_des_uniform_4_3060nonvar.HAL.290312/post_data';
% dirnameHFVar='../results/exp_des_uniform_4_3060.HAL.280312/post_data';
% 
% dirnameFullFact='../results/exp_des_3_3060.HAL.100212/data';

dirnameNatVar='../results/OHara2011_endo_uniform_30/post_data';
% dirnameHFNoVar='../results/OHara2011_endo_uniform_3060nonvar/post_data';
dirnameHFVar='../results/OHara2011_endo_uniform_3060/post_data';

dirnameFullFact='../results/OHara2011_endo_extremecase/data';

dirpost='../results/';
cd(dirpost);
%[sucess,message,messageid]=mkdir('statistics_correlations');
[sucess,message,messageid]=mkdir('statistics_correlations_Chaste');
if(sucess==0)
    fprintf('dir statistics_correlations. %s. Matlab error: %s\n',message,messageid);
    return;
elseif(isempty(message))
    fprintf('dir statistics_correlations. %s. Matlab error: %s\n',message,messageid);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

biomarker_string ={'APD_{80} (ms)' 'APD_{30}/APD_{80}' 'CaiTD_{80} (ms)' ...
    'CaiTD_{30}/CaiTD_{80}' 'CaiT_{max} (mMol)' 'DAPCaiT (ms)'};
%biomarker_string{2}

cd(curr_dir)
baseline=load([dirnameNatVar filesep 'baselineBio.dat']);
BCL = baseline(:,1);
APD80_base = baseline(:,2);
APD3080_base = baseline(:,3);
CaiTD80_base = baseline(:,4);
CaiTD3080_base = baseline(:,5);
DAPCaiT_base = baseline(:,6);
CaiTmax_base = baseline(:,7);
clear baselineBio;
nexp_baseline = 16385;
num_BCLs = length(BCL);
num_biomarkers = 6;

% orig = load([dirnameFullFact filesep 'DYN_2151.out']);
orig = load([dirnameFullFact filesep 'DYN_1.out']); % Chaste
APD80_60s   = orig(:,5)';
APD3080_60s = (orig(:,3)./orig(:,5))';
CaiTD80_60s = orig(:,9)';
CaiTD3080_60s = (orig(:,7)./orig(:,9))';
DAPCaiT_60s = orig(:,11)';
CaiTmax_60s = orig(:,13)';
%
% Generate matrix to hold all data:
allOutputsNatVar = zeros(4^7, num_BCLs, num_biomarkers);
% allOutputsHFNoVar = zeros(4^7, num_BCLs, num_biomarkers);
allOutputsHFVar = zeros(4^7, num_BCLs, num_biomarkers);

% Load explanatory data

HFExplanatoryData=load([dirnameHFVar filesep 'exp_design.dat']);
HFmask = HFExplanatoryData(1,:);
HFExplanatoryData = HFExplanatoryData(2:end-1,:);
shared_index_HF = get_intersection(HFExplanatoryData);

NVExplanatoryData=load([dirnameNatVar filesep 'exp_design.dat']);
NVmask = HFExplanatoryData(1,:);
NVExplanatoryData = NVExplanatoryData(2:end-1,:);
shared_index_NV = get_intersection(NVExplanatoryData);

%
% PROCESSING APD80
%
% APD80 data (Natural Variation)
%
A=load([dirnameNatVar filesep 'APD80.dat']);
% remove baseline point
A=A(1:4^7,:);
% load into data matrix:
allOutputsNatVar(:,:,1) = (A+ones(size(A))).*(ones(4^7,1)*APD80_base');
%
% APD80 data (Failing, No Natural Variation)
%
% 
% B=load([dirnameHFNoVar filesep 'APD80.dat']);
% % remove baseline point
% B=B(1:4^7,:);
% % load into data matrix:
% allOutputsHFNoVar(:,:,1) = (B+ones(size(B))).*(ones(4^7,1)*APD80_base');
% %
% APD80 data (Failing + Natural Variation)
%
C=load([dirnameHFVar filesep 'APD80.dat']);
% remove baseline point
C=C(1:4^7,:);
% load into data matrix:
allOutputsHFVar(:,:,1) = (C+ones(size(C))).*(ones(4^7,1)*APD80_base');
%
% PROCESSING APD3080
%
% APD3080 data (Natural Variation)
%

A=load([dirnameNatVar filesep 'APD3080.dat']);
% remove baseline point
A=A(1:4^7,:);
% load into data matrix:
allOutputsNatVar(:,:,2) = (A+ones(size(A))).*(ones(4^7,1)*APD3080_base');
%
% APD3080 data (Failing, No Natural Variation)
% %
% B=load([dirnameHFNoVar filesep 'APD3080.dat']);
% % remove baseline point
% B=B(1:4^7,:);
% % load into data matrix:
% allOutputsHFNoVar(:,:,2) = (B+ones(size(B))).*(ones(4^7,1)*APD3080_base');
% %
% APD3080 data (Failing + Natural Variation)
%
C=load([dirnameHFVar filesep 'APD3080.dat']);
% remove baseline point
C=C(1:4^7,:);
% load into data matrix:
allOutputsHFVar(:,:,2) = (C+ones(size(C))).*(ones(4^7,1)*APD3080_base');
%
% PROCESSING CaiTD80
%
% CaiTD80 data (Natural Variation)
%
A=load([dirnameNatVar filesep 'CaiTD80.dat']);
% remove baseline point
A=A(1:4^7,:);
% load into data matrix:
allOutputsNatVar(:,:,3) = (A+ones(size(A))).*(ones(4^7,1)*CaiTD80_base');
%
% CaiTD80 data (Failing, No Natural Variation)
% %
% B=load([dirnameHFNoVar filesep 'CaiTD80.dat']);
% % remove baseline point
% B=B(1:4^7,:);
% % load into data matrix:
% allOutputsHFNoVar(:,:,3) = (B+ones(size(B))).*(ones(4^7,1)*CaiTD80_base');
% %
% % CaiTD80 data (Failing + Natural Variation)
%
C=load([dirnameHFVar filesep 'CaiTD80.dat']);
% remove baseline point
C=C(1:4^7,:);
% load into data matrix:
allOutputsHFVar(:,:,3) = (C+ones(size(C))).*(ones(4^7,1)*CaiTD80_base');
%
% PROCESSING CaiRatio
%
% CaiTriangulation data (Hypothesis 1)
A=load([dirnameNatVar filesep 'Cai3080.dat']);
% remove baseline point
A=A(1:4^7,:);
% load into data matrix:
allOutputsNatVar(:,:,4) = (A+ones(size(A))).*(ones(4^7,1)*CaiTD3080_base');
%
% CaiTriangulation data (Failing, No Natural Variation)
% %
% B=load([dirnameHFNoVar filesep 'Cai3080.dat']);
% % remove baseline point
% B=B(1:4^7,:);
% % load into data matrix:
% allOutputsHFNoVar(:,:,4) = (B+ones(size(B))).*(ones(4^7,1)*CaiTD3080_base');
% %
% CaiTriangulation data (Failing + Natural Variation)
%
C=load([dirnameHFVar filesep 'Cai3080.dat']);
% remove baseline point
C=C(1:4^7,:);
% load into data matrix:
allOutputsHFVar(:,:,4) = (C+ones(size(C))).*(ones(4^7,1)*CaiTD3080_base');
%
% PROCESSING Systolic Cai
%
% Systolic Cai data (Natural Variation)
%
A=load([dirnameNatVar filesep 'SysCai.dat']);
% remove baseline point
A=A(1:4^7,:);
% load into data matrix:
allOutputsNatVar(:,:,5) = (A+ones(size(A))).*(ones(4^7,1)*CaiTmax_base')*1e4;
%
% % Systolic Cai data (Failing, No Natural Variation)
% %
% B=load([dirnameHFNoVar filesep 'SysCai.dat']);
% % remove baseline point
% B=B(1:4^7,:);
% % load into data matrix:
% allOutputsHFNoVar(:,:,5) = (B+ones(size(B))).*(ones(4^7,1)*CaiTmax_base')*1e4;
% %
% % Systolic Cai data (Failing + Natural Variation)
%
C=load([dirnameHFVar filesep 'SysCai.dat']);
% remove baseline point
C=C(1:4^7,:);
% load into data matrix:
allOutputsHFVar(:,:,5) = (C+ones(size(C))).*(ones(4^7,1)*CaiTmax_base')*1e4;
%
% PROCESSING AP - CaiT Delay
%
% AP-CaiT Delay data (Natural Variation)
%
A=load([dirnameNatVar filesep 'DAPCaiT.dat']);
% remove baseline point
A=A(1:4^7,:);
% load into data matrix:
allOutputsNatVar(:,:,6) = (A+ones(size(A))).*(ones(4^7,1)*DAPCaiT_base');
%
% AP-CaiT Delay data (Failing, No Natural Variation)
% %
% B=load([dirnameHFNoVar filesep 'DAPCaiT.dat']);
% % remove baseline point
% B=B(1:4^7,:);
% % load into data matrix:
% allOutputsHFNoVar(:,:,6) = (B+ones(size(B))).*(ones(4^7,1)*DAPCaiT_base');
% %
% AP-CaiT Delay data (Failing + Natural Variation)
%
C=load([dirnameHFVar filesep 'DAPCaiT.dat']);
% remove baseline point
C=C(1:4^7,:);
% load into data matrix:
allOutputsHFVar(:,:,6) = (C+ones(size(C))).*(ones(4^7,1)*DAPCaiT_base');
%


% Get Alternans locations

% load HF alternans data:
HF_APAlt = load([dirnameHFVar filesep 'APAlt.dat']);
[index_HFVar_APAlt_i,index_HFVar_APAlt_j]=find(HF_APAlt>0);
HF_CaiTAlt = load([dirnameHFVar filesep 'CaAlt.dat']);
[index_HFVar_CaAlt_i,index_HFVar_CaAlt_j]=find(HF_CaiTAlt>0);
% join AP and Ca Alternans

HFVar_alternans = union([index_HFVar_APAlt_i, index_HFVar_APAlt_j],[index_HFVar_CaAlt_i, index_HFVar_CaAlt_j], 'rows');

% % load HF no nat var alternans data
% HF_APAlt = load([dirnameHFNoVar filesep 'APAlt.dat']);
% [index_HFNoVar_APAlt_i,index_HFNoVar_APAlt_j]=find(HF_APAlt>0);
% HF_CaiTAlt = load([dirnameHFNoVar filesep 'CaAlt.dat']);
% [index_HFNoVar_CaAlt_i,index_HFNoVar_CaAlt_j]=find(HF_CaiTAlt>0);
% 
% HFNoVar_alternans = union([index_HFNoVar_APAlt_i, index_HFNoVar_APAlt_j],[index_HFNoVar_CaAlt_i, index_HFNoVar_CaAlt_j], 'rows');

% load NV alternans data
NV_APAlt = load([dirnameNatVar filesep 'APAlt.dat']);
[index_NV_APAlt_i,index_NV_APAlt_j]=find(NV_APAlt>0);
NV_CaiTAlt = load([dirnameNatVar filesep 'CaAlt.dat']);
[index_NV_CaAlt_i,index_NV_CaAlt_j]=find(NV_CaiTAlt>0);

NV_alternans = union([index_NV_APAlt_i, index_NV_APAlt_j],[index_NV_CaAlt_i, index_NV_CaAlt_j], 'rows');


if export_figures
    %cd([dirpost filesep 'statistics_correlations']);
    cd([dirpost filesep 'statistics_correlations_Chaste']);
end

PlotCorrelationPlots({allOutputsHFVar,allOutputsNatVar}, {HFVar_alternans,NV_alternans}, {HFExplanatoryData, NVExplanatoryData},BCL, biomarker_string, 'NonFailing and FailingVar', export_figures);
%PlotCorrelationPlots(allOutputsHFNoVar, BCL, biomarker_string, 'Failing', export_figures);
%PlotCorrelationPlots(allOutputsHFVar, BCL, biomarker_string, 'FailingVar', export_figures);

end

function PlotCorrelationPlots(matrices, alternans_indices, explanatory_data, BCL, biomarker_string, dataset_string, export_figures)
assert(length(matrices)==length(alternans_indices));
% Loop over cycle lengths
for cl = 1:length(BCL)
    
    fig=figure('Units','normalized','Position',[0.05 0.05 0.9 0.9]);
    % Loop over the different data sources
    for dataset_idx = 1:length(matrices)
        
        matrix = matrices{dataset_idx};
        explanatory = explanatory_data{dataset_idx};
        alternans = alternans_indices{dataset_idx};
        data_matrix = squeeze(matrix(:,cl,:));
        [explanatory_noalternans data_noalternans] = remove_alternans(explanatory,data_matrix,alternans,cl);
        shared_indices_noalternans = get_intersection(explanatory_noalternans);
        %size(data_noalternans)
        data{dataset_idx} = data_noalternans;
        shared{dataset_idx} = shared_indices_noalternans;
    end
    
    [H,AX,BigAx]=JW_plotmatrix_triangle(data{1},data{2}, shared{1}, shared{2});
    
    % Find the handles to all axes. Refer to Axes Properties Section of the
    % documentation for more detials on what each properties do.
    axes_handles = findobj(gcf,'Type','Axes');
    axes_posns = get(axes_handles(:),'Position');
    axes_posns_mat = cell2mat(axes_posns);
    xposns = sort(axes_posns_mat(:,1));
    yposns = sort(axes_posns_mat(:,2));
    xreq = [];
    yreq = [];
    
    % find the unique X-positions for the axes Position vector
    for i =1:length(xposns)
        p = find(xposns(i)==xposns(:));
        if length(p)>1 && xposns(i) ~= NaN
            xreq = [xreq;xposns(i)];
            xposns(p) = NaN;
        end
    end
    
    % find the unique Y-positions
    for i =1:length(yposns)
        p = find(yposns(i)==yposns(:));
        if length(p)>1 && yposns(i) ~= NaN
            yreq = [yreq;yposns(i)];
            yposns(p) = NaN;
        end
    end
    % Insert label or title
    for i = 1:length(axes_handles)
        
        % These make the background black...
        %background_colour = [0 0 0];
        %set(axes_handles(i),'Color',background_colour);
        
        p = get(axes_handles(i),'Position');
        if (p(1) == xreq(1)) && ~isempty(p(2) == yreq(:))
            set(gcf,'CurrentAxes',axes_handles(i));
            %i
            %biomarker_string{i}
%            if i~=1
                [TF, loc] = ismember(i, [1 2 4 7 11]);
                if TF
                ylabel(biomarker_string{loc+1});
                end
%            end
        end
        %if (p(2) == yreq(end)) && ~isempty(p(1) == xreq(:))
            set(gcf,'CurrentAxes',axes_handles(i));
%            if i~=1
                [TF, loc] = ismember(i, [11 12 13 14 15]);
                if TF
                xlabel(biomarker_string{loc});
                end
%            end
 %       end
%        if (p(1) == xreq(1)) && (p(2)== yreq(end))
%            q = get(get(axes_handles(i),'children'),'marker');
%             if q == '*'
%                 set(axes_handles(i),'visible','off')
%             end
%        end
    end
    
    %title_string=sprintf('Correlation between %s biomarkers,BCL %.0f ms', dataset_string,BCL(cl));
    %figTitle=title(BigAx,title_string);
    %titleLoc=get(figTitle, 'Position');
    %set(figTitle,'Position', titleLoc+[0 0.05 0]);
    
   
    hAll = findall(gcf);
    for idx = 1 : length(hAll)
        try
            set(hAll(idx),'fontsize',15);
        catch
            % nothing
        end
    end
    
    set(gcf, 'Color', 'w')
    
    if export_figures
        filename=sprintf('Correlation_%s_%0.f', dataset_string,BCL(cl));
        %export_fig(filename, '-eps', '-transparent');
        %export_fig(filename, '-png', '-transparent');
        %export_fig(filename, '-pdf', '-transparent');
        export_fig(filename, '-bmp', '-m3');
    end
    %pause

    close all
end

end


function [explanatory_noalternans, biomarker_noalternans] = remove_alternans(explanatory, biomarker, alternans_indices, current_CL_index)
index_alt_i = alternans_indices(:,1);
index_alt_j = alternans_indices(:,2);
explanatory_noalternans=explanatory;
biomarker_noalternans=biomarker;
if any(index_alt_j == current_CL_index)
    I=find(index_alt_j==current_CL_index);
    explanatory_noalternans(index_alt_i(I),:) = [];
    biomarker_noalternans(index_alt_i(I),:)=[];
end
end

% function indices_noalternans = remove_alternans_indices(indices,alternans,current_CL_index)
% index_alt_i = alternans(:,1);
% index_alt_j = alternans(:,2);
% indices_noalternans=indices;
% if any(index_alt_j == current_CL_index)
%     I=find(index_alt_j==current_CL_index);
%     alt_indices_current_CL = alternans(I);
%     for i=1:length(alt_indices_current_CL)
%         J = find(indices_noalternans==alt_indices_current_CL(i));
%         indices_noalternans(J) = [];
%     end
% end
% end

function indices = get_intersection(explanatory_data)
[num_experiments ~] = size(explanatory_data);
indices = [];
for idx=1:num_experiments
    if (explanatory_data(idx,1) >= 0.7 && explanatory_data(idx,1) <= 1.0 && ...
            explanatory_data(idx,2) >= 0.7 && explanatory_data(idx,2) <= 1.0 && ...
            explanatory_data(idx,3) >= 0.7 && explanatory_data(idx,3) <= 1.3 && ...
            explanatory_data(idx,4) >= 0.7 && explanatory_data(idx,4) <= 1.3 && ...
            explanatory_data(idx,5) >= 1.0 && explanatory_data(idx,5) <= 1.3 && ...
            explanatory_data(idx,6) >= 0.7 && explanatory_data(idx,6) <= 1.0 && ...
            explanatory_data(idx,7) >= 0.7 && explanatory_data(idx,7) <= 1.0)
        indices = [indices; idx];
    end
end
end

    