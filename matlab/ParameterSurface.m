function ParameterSurface

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
% dirnameHeartFailure='../results/exp_design_60_surf/post_data';
% dirnameNormalVar='../results/exp_design_30_surf/post_data';
dirnameHeartFailure='../results/OHara2011_endo_60_surf/post_data';
dirnameNormalVar='../results/OHara2011_endo_30_surf/post_data';
dirpost='../results/';
cd(dirpost);
% [sucess,message,messageid]=mkdir('parameter_surfaces');
[sucess,message,messageid]=mkdir('parameter_surfaces_Chaste');
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
HF_CaiT30Alt = load('Ca30Alt.dat');
[index_HF_Ca30Alt_i,index_HF_Ca30Alt_j]=find(HF_CaiT30Alt>0);
cd(curr_dir);
% load NV alternans data
cd(dirnameNormalVar);
NV_APAlt = load('APAlt.dat');
[index_NV_APAlt_i,index_NV_APAlt_j]=find(NV_APAlt>0);
NV_CaiTAlt = load('CaAlt.dat');
[index_NV_CaAlt_i,index_NV_CaAlt_j]=find(NV_CaiTAlt>0);
NV_CaiT30Alt = load('Ca30Alt.dat');
[index_NV_Ca30Alt_i,index_NV_Ca30Alt_j]=find(NV_CaiT30Alt>0);
cd(curr_dir)

% Now find the difference between calcium and AP alternans

HF_APAlt_indices = [index_HF_APAlt_i index_HF_APAlt_j];
HF_CaAlt_indices = [index_HF_CaAlt_i index_HF_CaAlt_j];
HF_Ca30Alt_indices = [index_HF_Ca30Alt_i index_HF_Ca30Alt_j];
NV_APAlt_indices = [index_NV_APAlt_i index_NV_APAlt_j];
NV_CaAlt_indices = [index_NV_CaAlt_i index_NV_CaAlt_j];
NV_Ca30Alt_indices = [index_NV_Ca30Alt_i index_NV_Ca30Alt_j];

cd(dirnameHeartFailure);
baseline=load('baselineBio.dat');
BCL = baseline(:,1);
CaiTD3080_baseline = baseline(:,5);
clear baselineBio;

% load HF parameter variations:
HFExplanatoryData=load('exp_design.dat');
mask = HFExplanatoryData(1,:);
HFExplanatoryData = HFExplanatoryData(2:end,:);

% load HF response variables
% Note, each column is for the different CLs contained in BCL
% APD80
HF_CaiTD3080 = load('Cai3080.dat');
clear Cai3080;

cd(curr_dir);
cd(dirnameNormalVar);

% load HF parameter variations:
NVExplanatoryData=load('exp_design.dat');
mask = NVExplanatoryData(1,:);
NVExplanatoryData = NVExplanatoryData(2:end,:);

% load HF response variables
% Note, each column is for the different CLs contained in BCL
% APD80
NV_CaiTD3080 = load('Cai3080.dat');
clear Cai3080;

cd(curr_dir);
cd(dirpost);
%cd('parameter_surfaces');
cd('parameter_surfaces_Chaste');

for i=1:length(BCL)
   
    % remove alternans points, use CaiTD80
    [NonFailingExplanatoryNoAlt,NonFailingCaiTD3080NoAlt] = remove_alternans(NVExplanatoryData, NV_CaiTD3080(:,i), NV_CaAlt_indices(:,1), NV_CaAlt_indices(:,2), i);
    [FailingExplanatoryNoAlt,FailingCaiTD3080NoAlt] = remove_alternans(HFExplanatoryData, HF_CaiTD3080(:,i), HF_CaAlt_indices(:,1), HF_CaAlt_indices(:,2), i);
 
    fig=figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    % fudge to ensure plot is on 3d axes
    plot3(NonFailingExplanatoryNoAlt(1,1),NonFailingExplanatoryNoAlt(1,2), (1+NonFailingCaiTD3080NoAlt(1))*baseline(i), 'k.');
    hold on
    for j=1:length(NonFailingCaiTD3080NoAlt)
        if (NonFailingExplanatoryNoAlt(j,1) <= 1.0 && NonFailingExplanatoryNoAlt(j,1) >= 0.7 ...
                && NonFailingExplanatoryNoAlt(j,2) <= 1.0 && NonFailingExplanatoryNoAlt(j,2) >= 0.7)
            plot3(NonFailingExplanatoryNoAlt(j,1),NonFailingExplanatoryNoAlt(j,2), (1+NonFailingCaiTD3080NoAlt(j))*CaiTD3080_baseline(i), 'color', 0.2*[1 1 1], 'marker', '.', 'MarkerSize', 30);
        else 
            plot3(NonFailingExplanatoryNoAlt(j,1),NonFailingExplanatoryNoAlt(j,2), (1+NonFailingCaiTD3080NoAlt(j))*CaiTD3080_baseline(i), 'color', 0.45*[1 1 1], 'marker', '.', 'MarkerSize', 30);
        end
    end
    for j=1:length(FailingCaiTD3080NoAlt)
        if (FailingExplanatoryNoAlt(j,1) <= 1.0 && FailingExplanatoryNoAlt(j,1) >= 0.7 ...
                && FailingExplanatoryNoAlt(j,2) <= 1.0 && FailingExplanatoryNoAlt(j,2) >= 0.7)
            plot3(FailingExplanatoryNoAlt(j,1),FailingExplanatoryNoAlt(j,2), (1+FailingCaiTD3080NoAlt(j))*CaiTD3080_baseline(i), 'color', 0.2*[1 1 1], 'marker', '.', 'MarkerSize', 30);
        else 
            plot3(FailingExplanatoryNoAlt(j,1),FailingExplanatoryNoAlt(j,2), (1+FailingCaiTD3080NoAlt(j))*CaiTD3080_baseline(i), 'color', 0.7*[1 1 1], 'marker', '.', 'MarkerSize', 30);
        end
    end
    hold off
    view([50, 24]);
    set(gca,'zlim', [0.35 0.65]);
    
    xlabel(mask_string{mask(1)+1});
    ylabel(mask_string{mask(2)+1});
    zlabel('CaiTD3080');
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
    
    % save to file using export_fig
    sname=sprintf('ParameterSurface_%.0f', BCL(i));
    export_fig([sname '.png'], '-png', '-transparent');
    export_fig([sname '.eps'], '-eps', '-transparent');
    close all

    
end

cd(curr_dir);

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

function normalised_data = Normalise_Data(data)
[data_rows, ~] = size(data);
normalised_data = (data - (ones(data_rows,1)*mean(data)))./(ones(data_rows,1)*sqrt(var(data)));
end