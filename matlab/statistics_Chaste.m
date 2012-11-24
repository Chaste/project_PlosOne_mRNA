%
% Computes the statistics and histograms for the following biomarkers:
% APD80; AP triangulation as ADP80-ADP30; CaiTD;
% Cai triangulation as CaiTD80-CaiTD30; AP Cai delay
%
% Data is normalized with respect to the baseline value
%
% postprocessing program
%
function statistics_boxplot
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
% dirnameNatVar='../results/exp_des_uniform_4_30.HAL.260312/post_data';
% dirnameHFNoVar='../results/exp_des_uniform_4_3060nonvar.HAL.290312/post_data';
% dirnameHFVar='../results/exp_des_uniform_4_3060.HAL.280312/post_data';
% 
% dirnameFullFact='../results/exp_des_3_3060.HAL.100212/data';

dirnameNatVar='../results/OHara2011_endo_uniform_30/post_data';
%dirnameHFNoVar='../results/OHara2011_endo_uniform_3060nonvar/post_data';
dirnameHFVar='../results/OHara2011_endo_uniform_3060/post_data';

dirnameFullFact='../results/OHara2011_endo_extremecase/data';


dirpost='../results/';
cd(dirpost);
%[sucess,message,messageid]=mkdir('statistics_boxplot');
[sucess,message,messageid]=mkdir('statistics_boxplot_Chaste');
% if(sucess==0)
%     fprintf('dir statistics_boxplot. %s. Matlab error: %s\n',message,messageid);
%     return;
% elseif(isempty(message))
%     fprintf('dir statistics_boxplot. %s. Matlab error: %s\n',message,messageid);
% end
if(sucess==0)
    fprintf('dir statistics_boxplot_Chaste. %s. Matlab error: %s\n',message,messageid);
    return;
elseif(isempty(message))
    fprintf('dir statistics_boxplot_Chaste. %s. Matlab error: %s\n',message,messageid);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
cd(curr_dir)
cd(dirnameNatVar);
baseline=load('baselineBio.dat');
BCL = baseline(:,1);
APD80_base = baseline(:,2);
APD3080_base = baseline(:,3);
CaiTD80_base = baseline(:,4);
CaiTD3080_base = baseline(:,5);
DAPCaiT_base = baseline(:,6);
CaiTmax_base = baseline(:,7);
CaiTD30_base = baseline(:,11);
clear baselineBio;
nexp_baseline = 16385;

cd(curr_dir)
cd(dirnameFullFact);

% orig = load('DYN_2151.out'); 
orig = load('DYN_1.out'); % Change needed for Chaste
APD80_60s   = orig(:,5)';
APD3080_60s = (orig(:,3)./orig(:,5))';
for i=1:length(BCL)
    if APD80_60s(i)==-1
        APD3080_60s=-1;
    end
end
CaiTD80_60s = orig(:,9)';
CaiTD3080_60s = (orig(:,7)./orig(:,9))';
for i=1:length(BCL)
    if CaiTD80_60s(i)==-1
        CaiTD3080_60s=-1;
    end
end
DAPCaiT_60s = orig(:,11)';
CaiTmax_60s = orig(:,13)';
CaiTD30_60s = orig(:,7)';
cd(curr_dir)
%
% loading experimental design. HF
%
% ffH1=load('exp_design.dat');
% mask = ffH1(1,:);
% ffH1 = ffH1(2:end,:);
% nexp_baseH1 = find(prod(ffH1,2)==1);
% nexp=length(ffH1(:,1));
% Id = find((ffH1(:,3)==1)&(ffH1(:,4)==1));
% i1=find(Id==nexp_baseH1);
% Id(i1)=0;
% i1=find(Id>0);
% IdH1=Id(i1);
% IdH2=[1:(nexp_baseH1-1) (nexp_baseH1+1):nexp];
%
%

% Get Alternans locations
% load HF alternans data:
cd(dirnameHFVar)
HF_APAlt = load('APAlt.dat');
[index_HFVar_APAlt_i,index_HFVar_APAlt_j]=find(HF_APAlt>0);
HF_CaiTAlt = load('CaAlt.dat');
[index_HFVar_CaAlt_i,index_HFVar_CaAlt_j]=find(HF_CaiTAlt>0);
HF_CaiT30Alt = load('Ca30Alt.dat');
[index_HFVar_Ca30Alt_i,index_HFVar_Ca30Alt_j]=find(HF_CaiT30Alt>0);
cd(curr_dir);
% load HF no nat var alternans data
% cd(dirnameHFNoVar)
% HF_APAlt = load('APAlt.dat');
% [index_HFNoVar_APAlt_i,index_HFNoVar_APAlt_j]=find(HF_APAlt>0);
% HF_CaiTAlt = load('CaAlt.dat');
% [index_HFNoVar_CaAlt_i,index_HFNoVar_CaAlt_j]=find(HF_CaiTAlt>0);
% HF_CaiT30Alt = load('Ca30Alt.dat');
% [index_HFNoVar_Ca30Alt_i,index_HFNoVar_Ca30Alt_j]=find(HF_CaiT30Alt>0);
% cd(curr_dir);
% load NV alternans data
cd(dirnameNatVar);
NV_APAlt = load('APAlt.dat');
[index_NV_APAlt_i,index_NV_APAlt_j]=find(NV_APAlt>0);
NV_CaiTAlt = load('CaAlt.dat');
[index_NV_CaAlt_i,index_NV_CaAlt_j]=find(NV_CaiTAlt>0);
NV_CaiT30Alt = load('CaAlt.dat');
[index_NV_Ca30Alt_i,index_NV_Ca30Alt_j]=find(NV_CaiT30Alt>0);
cd(curr_dir)

HFVar_APAlt_indices = [index_HFVar_APAlt_i index_HFVar_APAlt_j];
HFVar_CaAlt_indices = [index_HFVar_CaAlt_i index_HFVar_CaAlt_j];
HFVar_Ca30Alt_indices = [index_HFVar_Ca30Alt_i index_HFVar_Ca30Alt_j];
% HFNoVar_APAlt_indices = [index_HFNoVar_APAlt_i index_HFNoVar_APAlt_j];
% HFNoVar_CaAlt_indices = [index_HFNoVar_CaAlt_i index_HFNoVar_CaAlt_j];
% HFNoVar_Ca30Alt_indices = [index_HFNoVar_Ca30Alt_i index_HFNoVar_Ca30Alt_j];
NV_APAlt_indices = [index_NV_APAlt_i index_NV_APAlt_j];
NV_CaAlt_indices = [index_NV_CaAlt_i index_NV_CaAlt_j];
NV_Ca30Alt_indices = [index_NV_Ca30Alt_i index_NV_Ca30Alt_j];
% PROCESSING APD80
%
% APD80 data (Natural Variation)
%
cd(curr_dir);
cd(dirnameNatVar);
A=load('APD80.dat');
% make into percentages
A=100*A;
%
% APD80 data (Failing, No Natural Variation)
%
% cd(curr_dir);
% cd(dirnameHFNoVar);
% B=load('APD80.dat');
% % make into percentages
% B=100*B;
%
% APD80 data (Failing + Natural Variation)
%
cd(curr_dir);
cd(dirnameHFVar);
C=load('APD80.dat');
% make into percentages
C=100*C;
%
cd(curr_dir);
cd(dirpost);
% cd('statistics_boxplot');
cd('statistics_boxplot_Chaste');
%
% Hypothesis 1
%
% fp=fopen('APD80H1.dat','w');
% fprintf(fp,'%% values represent represent relative change with respect to\n');
% fprintf(fp,'%% control value given by OHara base model\n');
% fprintf(fp,'%%1st row BCL\n');
% fprintf(fp,'%%2nd row baseline\n');
% fprintf(fp,'%%3rd row average without baseline point\n');
% fprintf(fp,'%%4th row standard deviation without baseline point\n');
% fprintf(fp,'%6.1f ',BCL);
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',B(nexp_baseline,:));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',mean(B(1:end-1,:)));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',std(B(1:end-1,:)));
% fprintf(fp,'\n');
% fclose(fp);
%
% Hypothesis 2
%
fp=fopen('APD80H2.dat','w');
fprintf(fp,'%% values represent represent relative change with respect to\n');
fprintf(fp,'%% control value given by OHara base model\n');
fprintf(fp,'%%1st row BCL\n');
fprintf(fp,'%%2nd row average physiologic\n');
fprintf(fp,'%%2nd row standard deviation physiologic\n');
fprintf(fp,'%%3rd row average without baseline point\n');
fprintf(fp,'%%4th row standard deviation without baseline point\n');
fprintf(fp,'%6.1f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(C(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(C(1:end-1,:)));
fprintf(fp,'\n');
fclose(fp);

% boxplots
%PlotBoxPlots(BCL, A,B,C, APD80_base, APD80_60s,NV_APAlt_indices,HFNoVar_APAlt_indices,HFVar_APAlt_indices,'APD80', 'APD_{80}','(ms)')
% histograms
% PlotHistograms(BCL, A,C, APD80_base, APD80_60s,NV_APAlt_indices,HFVar_APAlt_indices,'APD80', 'APD_{80}','(ms)')
% APD80limits = [110 490];
% PlotHistogramsActualValuesOnly(BCL, A,C, APD80_base, APD80_60s,NV_APAlt_indices,HFVar_APAlt_indices,'APD80', 'APD_{80}','(ms)')

%
%
% PROCESSING APDTriangulation
%
% APDTriangulation data (Hypothesis 1)
%
%
% APD3080 data (Natural Variation)
%
cd(curr_dir);
cd(dirnameNatVar);
A=load('APD3080.dat');
% make into percentages
A=100*A;
%
% APD3080 data (Failing, No Natural Variation)
% %
% cd(curr_dir);
% cd(dirnameHFNoVar);
% B=load('APD3080.dat');
% % make into percentages
% B=100*B;
%
% APD3080 data (Failing + Natural Variation)
%
cd(curr_dir);
cd(dirnameHFVar);
C=load('APD3080.dat');
% make into percentages
C=100*C;
%
cd(curr_dir);
cd(dirpost);
%cd('statistics_boxplot');
cd('statistics_boxplot_Chaste');
%
% Hypothesis 1
%
% fp=fopen('APD3080H1.dat','w');
% fprintf(fp,'%% values represent represent relative change with respect to\n');
% fprintf(fp,'%% control value given by OHara base model\n');
% fprintf(fp,'%%1st row BCL\n');
% fprintf(fp,'%%2nd row baseline\n');
% fprintf(fp,'%%3rd row average without baseline point\n');
% fprintf(fp,'%%4th row standard deviation without baseline point\n');
% fprintf(fp,'%6.1f ',BCL);
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',B(nexp_baseline,:));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',mean(B(1:end-1,:)));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',std(B(1:end-1,:)));
% fprintf(fp,'\n');
% fclose(fp);
% %
% Hypothesis 2
%
fp=fopen('APD3080H2.dat','w');
fprintf(fp,'%% values represent represent relative change with respect to\n');
fprintf(fp,'%% control value given by OHara base model\n');
fprintf(fp,'%%1st row BCL\n');
fprintf(fp,'%%2nd row average physiologic\n');
fprintf(fp,'%%2nd row standard deviation physiologic\n');
fprintf(fp,'%%3rd row average without baseline point\n');
fprintf(fp,'%%4th row standard deviation without baseline point\n');
fprintf(fp,'%6.1f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(C(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(C(1:end-1,:)));
fprintf(fp,'\n');
fclose(fp);

% boxplots
%PlotBoxPlots(BCL, A,B,C, APD3080_base, APD3080_60s,NV_APAlt_indices,HFNoVar_APAlt_indices,HFVar_APAlt_indices,'APD3080', 'APD_{30}/APD_{80}','')
% histograms
% PlotHistograms(BCL, A,C, APD3080_base, APD3080_60s,NV_APAlt_indices,HFVar_APAlt_indices,'APD3080', 'APD_{30}/APD_{80}','')
% APD3080limits=[0.48 0.72];
% PlotHistogramsActualValuesOnly(BCL, A,C, APD3080_base, APD3080_60s,NV_APAlt_indices,HFVar_APAlt_indices,'APD3080', 'APD_{30}/APD_{80}','')
%
%
% PROCESSING CaiTD80
%
% CaiTD80 data (Natural Variation)
%
cd(curr_dir);
cd(dirnameNatVar);
A=load('CaiTD80.dat');
% make into percentages
A=100*A;
%
% CaiTD80 data (Failing, No Natural Variation)
%
% cd(curr_dir);
% cd(dirnameHFNoVar);
% B=load('CaiTD80.dat');
% % make into percentages
% B=100*B;
% %
% CaiTD80 data (Failing + Natural Variation)
%
cd(curr_dir);
cd(dirnameHFVar);
C=load('CaiTD80.dat');
% make into percentages
C=100*C;
%
cd(curr_dir);
cd(dirpost);
% cd('statistics_boxplot');
cd('statistics_boxplot_Chaste');
%
% Hypothesis 1
%
% fp=fopen('CaiTD80H1.dat','w');
% fprintf(fp,'%% values represent represent relative change with respect to\n');
% fprintf(fp,'%% control value given by OHara base model\n');
% fprintf(fp,'%%1st row BCL\n');
% fprintf(fp,'%%2nd row baseline\n');
% fprintf(fp,'%%3rd row average without baseline point\n');
% fprintf(fp,'%%4th row standard deviation without baseline point\n');
% fprintf(fp,'%6.1f ',BCL);
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',B(nexp_baseline,:));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',mean(B(1:end-1,:)));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',std(B(1:end-1,:)));
% fprintf(fp,'\n');
% fclose(fp);
%
% Hypothesis 2
%
fp=fopen('CaiTD80H2.dat','w');
fprintf(fp,'%% values represent represent relative change with respect to\n');
fprintf(fp,'%% control value given by OHara base model\n');
fprintf(fp,'%%1st row BCL\n');
fprintf(fp,'%%2nd row average physiologic\n');
fprintf(fp,'%%2nd row standard deviation physiologic\n');
fprintf(fp,'%%3rd row average without baseline point\n');
fprintf(fp,'%%4th row standard deviation without baseline point\n');
fprintf(fp,'%6.1f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(C(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(C(1:end-1,:)));
fprintf(fp,'\n');
fclose(fp);
CaiTD80_60s
% boxplots
%PlotBoxPlots(BCL, A,B,C, CaiTD80_base, CaiTD80_60s,NV_CaAlt_indices,HFNoVar_CaAlt_indices,HFVar_CaAlt_indices,'CaiTD80', 'CaiTD_{80}','(ms)')
% histograms
% PlotHistograms(BCL, A,C, CaiTD80_base, CaiTD80_60s, NV_CaAlt_indices,HFVar_CaAlt_indices,'CaiTD80', 'CaiTD_{80}','(ms)')
% CaiTD80limits = [160 990];
% PlotHistogramsActualValuesOnly(BCL, A,C, CaiTD80_base, CaiTD80_60s, NV_CaAlt_indices,HFVar_CaAlt_indices,'CaiTD80', 'CaiTD_{80}','(ms)')
%

%
%
% PROCESSING CaiTD30
%
% CaiTD30 data (Natural Variation)
%
cd(curr_dir);
cd(dirnameNatVar);
A=load('CaiTD30.dat');
% make into percentages
A=100*A;
%
% CaiTD30 data (Failing, No Natural Variation)
%
% cd(curr_dir);
% cd(dirnameHFNoVar);
% B=load('CaiTD30.dat');
% % make into percentages
% B=100*B;
% %
% CaiTD30 data (Failing + Natural Variation)
%
cd(curr_dir);
cd(dirnameHFVar);
C=load('CaiTD30.dat');
% make into percentages
C=100*C;
%
cd(curr_dir);
cd(dirpost);
%cd('statistics_boxplot');
cd('statistics_boxplot_Chaste');
%
% % Hypothesis 1
% %
% fp=fopen('CaiTD30H1.dat','w');
% fprintf(fp,'%% values represent represent relative change with respect to\n');
% fprintf(fp,'%% control value given by OHara base model\n');
% fprintf(fp,'%%1st row BCL\n');
% fprintf(fp,'%%2nd row baseline\n');
% fprintf(fp,'%%3rd row average without baseline point\n');
% fprintf(fp,'%%4th row standard deviation without baseline point\n');
% fprintf(fp,'%6.1f ',BCL);
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',B(nexp_baseline,:));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',mean(B(1:end-1,:)));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',std(B(1:end-1,:)));
% fprintf(fp,'\n');
% fclose(fp);
% %
% Hypothesis 2
%
fp=fopen('CaiTD30H2.dat','w');
fprintf(fp,'%% values represent represent relative change with respect to\n');
fprintf(fp,'%% control value given by OHara base model\n');
fprintf(fp,'%%1st row BCL\n');
fprintf(fp,'%%2nd row average physiologic\n');
fprintf(fp,'%%2nd row standard deviation physiologic\n');
fprintf(fp,'%%3rd row average without baseline point\n');
fprintf(fp,'%%4th row standard deviation without baseline point\n');
fprintf(fp,'%6.1f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(C(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(C(1:end-1,:)));
fprintf(fp,'\n');
fclose(fp);

% boxplots
%PlotBoxPlots(BCL, A,B,C, CaiTD30_base, CaiTD30_60s,NV_CaAlt_indices,HFNoVar_CaAlt_indices,HFVar_CaAlt_indices,'CaiTD30', 'CaiTD_{80}','(ms)')
% histograms
% PlotHistograms(BCL, A,C, CaiTD30_base, CaiTD30_60s, NV_CaAlt_indices,HFVar_CaAlt_indices,'CaiTD30', 'CaiTD_{80}','(ms)')
%CaiTD30limits=[40 440];
% PlotHistogramsActualValuesOnly(BCL, A,C, CaiTD30_base, CaiTD30_60s, NV_Ca30Alt_indices,HFVar_Ca30Alt_indices,'CaiTD30', 'CaiTD_{30}','(ms)')
%



%
% PROCESSING CaiRatio
%
% CaiTriangulation data (Hypothesis 1)
cd(curr_dir);
cd(dirnameNatVar);
A=load('Cai3080.dat');
% make into percentages
A=100*A;
%
% CaiTriangulation data (Failing, No Natural Variation)
%
% cd(curr_dir);
% cd(dirnameHFNoVar);
% B=load('Cai3080.dat');
% % make into percentages
% B=100*B;
%
% CaiTriangulation data (Failing + Natural Variation)
%
cd(curr_dir);
cd(dirnameHFVar);
C=load('Cai3080.dat');
% make into percentages
C=100*C;
%
cd(curr_dir);
cd(dirpost);
% cd('statistics_boxplot');
cd('statistics_boxplot_Chaste');
%
% Hypothesis 1
%
% fp=fopen('Cai3080H1.dat','w');
% fprintf(fp,'%% values represent represent relative change with respect to\n');
% fprintf(fp,'%% control value given by OHara base model\n');
% fprintf(fp,'%%1st row BCL\n');
% fprintf(fp,'%%2nd row baseline\n');
% fprintf(fp,'%%3rd row average without baseline point\n');
% fprintf(fp,'%%4th row standard deviation without baseline point\n');
% fprintf(fp,'%6.1f ',BCL);
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',B(nexp_baseline,:));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',mean(B(1:end-1,:)));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',std(B(1:end-1,:)));
% fprintf(fp,'\n');
% fclose(fp);
% %
% Hypothesis 2
%
fp=fopen('Cai3080H2.dat','w');
fprintf(fp,'%% values represent represent relative change with respect to\n');
fprintf(fp,'%% control value given by OHara base model\n');
fprintf(fp,'%%1st row BCL\n');
fprintf(fp,'%%2nd row average physiologic\n');
fprintf(fp,'%%2nd row standard deviation physiologic\n');
fprintf(fp,'%%3rd row average without baseline point\n');
fprintf(fp,'%%4th row standard deviation without baseline point\n');
fprintf(fp,'%6.1f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(C(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(C(1:end-1,:)));
fprintf(fp,'\n');
fclose(fp);

% boxplots
% PlotBoxPlots(BCL, A,B,C, CaiTD3080_base, CaiTD3080_60s,NV_CaAlt_indices,HFNoVar_CaAlt_indices,HFVar_CaAlt_indices, 'CaiTD3080', 'CaiTD_{30}/CaiTD_{80}','')
% histograms
% PlotHistograms(BCL, A,C, CaiTD3080_base, CaiTD3080_60s,NV_CaAlt_indices,HFVar_CaAlt_indices, 'CaiTD3080', 'CaiTD_{30}/CaiTD_{80}','')
% CaiTD3080limits=[0 0.73];
% PlotHistogramsActualValuesOnly(BCL, A,C, CaiTD3080_base, CaiTD3080_60s,NV_CaAlt_indices,HFVar_CaAlt_indices, 'CaiTD3080', 'CaiTD_{30}/CaiTD_{80}','')

%
%
% PROCESSING Systolic Cai
%
% Systolic Cai data (Natural Variation)
%
cd(curr_dir);
cd(dirnameNatVar);
A=load('SysCai.dat');
% make into percentages
A=100*A;
%
% Systolic Cai data (Failing, No Natural Variation)
% %
% cd(curr_dir);
% cd(dirnameHFNoVar);
% B=load('SysCai.dat');
% % make into percentages
% B=100*B;
%
% Systolic Cai data (Failing + Natural Variation)
%
cd(curr_dir);
cd(dirnameHFVar);
C=load('SysCai.dat');
% make into percentages
C=100*C;
%
cd(curr_dir);
cd(dirpost);
% cd('statistics_boxplot');
cd('statistics_boxplot_Chaste');
%
% Hypothesis 1
%
% fp=fopen('SysCaiH1.dat','w');
% fprintf(fp,'%% values represent represent relative change with respect to\n');
% fprintf(fp,'%% control value given by OHara base model\n');
% fprintf(fp,'%%1st row BCL\n');
% fprintf(fp,'%%2nd row baseline\n');
% fprintf(fp,'%%3rd row average without baseline point\n');
% fprintf(fp,'%%4th row standard deviation without baseline point\n');
% fprintf(fp,'%6.1f ',BCL);
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',B(nexp_baseline,:));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',mean(B(1:end-1,:)));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',std(B(1:end-1,:)));
% fprintf(fp,'\n');
% fclose(fp);
%
% Hypothesis 2
%
fp=fopen('SysCaiH2.dat','w');
fprintf(fp,'%% values represent represent relative change with respect to\n');
fprintf(fp,'%% control value given by OHara base model\n');
fprintf(fp,'%%1st row BCL\n');
fprintf(fp,'%%2nd row average physiologic\n');
fprintf(fp,'%%2nd row standard deviation physiologic\n');
fprintf(fp,'%%3rd row average without baseline point\n');
fprintf(fp,'%%4th row standard deviation without baseline point\n');
fprintf(fp,'%6.1f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(C(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(C(1:end-1,:)));
fprintf(fp,'\n');
fclose(fp);

% boxplots
% PlotBoxPlots(BCL, A,B,C, CaiTmax_base, CaiTmax_60s,NV_CaAlt_indices,HFNoVar_CaAlt_indices,HFVar_CaAlt_indices, 'CaiTMax', 'CaiT_{max}', '(mMol)')
% histograms
% PlotHistograms(BCL, A,C, CaiTmax_base, CaiTmax_60s,NV_CaAlt_indices,HFVar_CaAlt_indices, 'CaiTMax', 'CaiT_{max}','(mMol)')
% CaiTMaxlimits = [0 0.0012];
% PlotHistogramsActualValuesOnly(BCL, A,C, CaiTmax_base, CaiTmax_60s,NV_CaAlt_indices,HFVar_CaAlt_indices, 'CaiTMax', 'CaiT_{max}','(mMol)')
%
%
% PROCESSING AP - CaiT Delay
%
% AP-CaiT Delay data (Natural Variation)
%
cd(curr_dir);
cd(dirnameNatVar);
A=load('DAPCaiT.dat');
% make into percentages
A=100*A;
% %
% % AP-CaiT Delay data (Failing, No Natural Variation)
% %
% cd(curr_dir);
% cd(dirnameHFNoVar);
% B=load('DAPCaiT.dat');
% % make into percentages
% B=100*B;
%
% AP-CaiT Delay data (Failing + Natural Variation)
%
cd(curr_dir);
cd(dirnameHFVar);
C=load('DAPCaiT.dat');
% make into percentages
C=100*C;
%
cd(curr_dir);
cd(dirpost);
% cd('statistics_boxplot');
cd('statistics_boxplot_Chaste');
%
% Hypothesis 1
%
% fp=fopen('DAPCaiTH1.dat','w');
% fprintf(fp,'%% values represent represent relative change with respect to\n');
% fprintf(fp,'%% control value given by OHara base model\n');
% fprintf(fp,'%%1st row BCL\n');
% fprintf(fp,'%%2nd row baseline\n');
% fprintf(fp,'%%3rd row average without baseline point\n');
% fprintf(fp,'%%4th row standard deviation without baseline point\n');
% fprintf(fp,'%6.1f ',BCL);
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',B(nexp_baseline,:));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',mean(B(1:end-1,:)));
% fprintf(fp,'\n');
% fprintf(fp,'%12.4e ',std(B(1:end-1,:)));
% fprintf(fp,'\n');
% fclose(fp);
%
% Hypothesis 2
%
fp=fopen('DAPCaiTH2.dat','w');
fprintf(fp,'%% values represent represent relative change with respect to\n');
fprintf(fp,'%% control value given by OHara base model\n');
fprintf(fp,'%%1st row BCL\n');
fprintf(fp,'%%2nd row average physiologic\n');
fprintf(fp,'%%2nd row standard deviation physiologic\n');
fprintf(fp,'%%3rd row average without baseline point\n');
fprintf(fp,'%%4th row standard deviation without baseline point\n');
fprintf(fp,'%6.1f ',BCL);
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(A(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',mean(C(1:end-1,:)));
fprintf(fp,'\n');
fprintf(fp,'%12.4e ',std(C(1:end-1,:)));
fprintf(fp,'\n');
fclose(fp);

% boxplots
%PlotBoxPlots(BCL, A,B,C, DAPCaiT_base, DAPCaiT_60s,NV_APAlt_indices,HFNoVar_APAlt_indices,HFVar_APAlt_indices, 'DAPCaiT', 'DAPCaiT','(ms)')
% histograms
% PlotHistograms(BCL, A,C, DAPCaiT_base, DAPCaiT_60s,NV_APAlt_indices,HFVar_APAlt_indices, 'DAPCaiT', 'DAPCaiT','(ms)')
% DAPCaiTlimits = [0 40];
PlotHistogramsActualValuesOnly(BCL, A,C, DAPCaiT_base, DAPCaiT_60s,NV_APAlt_indices,HFVar_APAlt_indices, 'DAPCaiT', 'DAPCaiT','(ms)')

cd(curr_dir);

end

function PlotBoxPlots(BCL, NonFailing, Failing, FailingWithVar, baseline, results_60s, NonFailingAltIndices, FailingAltIndices, FailingWithVarAltIndices, FileNameString, FigTitleString,units)
% note: don't need explanatory data values for this!
exp_data=zeros(16385,11);
for i=1:length(BCL)
    fig=figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    
    [~,NonFailingNoAlt] = remove_alternans(exp_data, NonFailing(:,i), NonFailingAltIndices(:,1), NonFailingAltIndices(:,2), i);
    [~,FailingNoAlt] = remove_alternans(exp_data, Failing(:,i), FailingAltIndices(:,1), FailingAltIndices(:,2), i);
    [~,FailingWithVarNoAlt] = remove_alternans(exp_data, FailingWithVar(:,i), FailingWithVarAltIndices(:,1), FailingWithVarAltIndices(:,2), i);
    % set groups
    % group = [repmat({'Non-Failing'}, size(B(IdH2, i))); repmat({'Failing'}, size(A(IdH1,i))); repmat({'Failing + Variability'}, size(A(IdH2,i)))];
    group = [repmat({'Non-Failing'}, size(NonFailingNoAlt)); repmat({'Failing'}, size(FailingNoAlt)); repmat({'Failing + Variability'}, size(FailingWithVarNoAlt))];
    % make box plots
    %bhandle=boxplot([B(IdH2, i); A(IdH1,i); A(IdH2,i)],group,'color','k');
    bhandle=boxplot([NonFailingNoAlt; FailingNoAlt; FailingWithVarNoAlt],group,'color','k');
    ylabel('% Change From Baseline');
    title_text=sprintf('%s, BCL=%.0f ms',FigTitleString,BCL(i));
    title(title_text);

    % Set line width of box plots
    set(bhandle,'linewidth',2.5);
        
    ax1 = gca;
    
    % loop over box plots, setting background colour
    h = findobj(ax1,'Tag','Box');
    patch(get(h(3),'XData'),get(h(3),'YData'),'b','FaceAlpha',.5);
    patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.3);
    patch(get(h(1),'XData'),get(h(1),'YData'),'r','FaceAlpha',.7);
    
    
    % Now add in second axis on the right:
    ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
    % Set units for the new axis using baseline values
    axis_limits = ylim(ax1);
    axis_limits_factor =[1 1]+axis_limits/100;
    new_limits =axis_limits_factor*baseline(i);
    ylim(ax2, new_limits);
    xlim(ax2, [0 1]);
    % Set width of axes
    set(ax1, 'linewidth', 2);
    set(ax2, 'linewidth', 2);
    
    % make background transparent
    set(ax1,'color','none');
    set(ax2,'color','none');
    
    % remove x axis
    set(ax1,'xcolor',get(fig,'color'),'xtick',[]);
    set(ax2,'xcolor',get(fig,'color'),'xtick',[]);
    
    % remove box
    set(ax1,'Box','off');
    set(ax2,'Box','off');
    
    % plot values for 60s case
    hold on
    plot(ax2, 0.166, results_60s(i), 'kx', 0.5, results_60s(i),'kx', 0.834, results_60s(i),'kx',...
        'MarkerSize', 20, 'LineWidth', 4);
    hold off
    
    % label second axis
    actual_value_label = sprintf('%s %s',FigTitleString, units);
    ylabel(ax2, actual_value_label);
    ylabh = get(ax2,'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') + [.05 0 0]);
    
    % set font size for all text
    hAll = findall(gcf);
    for idx = 1 : length(hAll)
        try
            set(hAll(idx),'fontsize',30);
        catch
            % nothing
        end
    end
       
    sname=sprintf('boxplot_%s_BCL%.0f',FileNameString, BCL(i));
    export_fig([sname '.png'], '-png', '-transparent');
    %print('-depsc',sname);
    close all
end
end

function PlotHistograms(BCL, NonFailing, FailingWithVar, baseline, results_60s, NonFailingAltIndices, FailingWithVarAltIndices, FileNameString, FigTitleString, units)
% note: don't need explanatory data values for this!
exp_data=zeros(16385,11);
for i=1:length(BCL)
    [~,NonFailingNoAlt] = remove_alternans(exp_data, NonFailing(:,i), NonFailingAltIndices(:,1), NonFailingAltIndices(:,2), i);
    [~,FailingWithVarNoAlt] = remove_alternans(exp_data, FailingWithVar(:,i), FailingWithVarAltIndices(:,1), FailingWithVarAltIndices(:,2), i);
    nbins = fix(sqrt(length(NonFailingNoAlt)));
    fig=figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
    hold on
    hist(FailingWithVarNoAlt,nbins);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75);
    hist(NonFailingNoAlt,nbins);
    h = findobj(gca,'Type','patch');
    set(h,'facealpha',0.75);
    hold off
    ax2 = gca;
    
    %current_position=get(ax2,'Position');
    %set(ax2, 'Position',current_position+[0 -0.015 0 -0.1]);
    
    histleg=legend('HF','Normal');
    set(histleg, 'Position', get(histleg, 'Position')+[-0.025 -0.1 0 0]);
    xlabel('% Change from baseline','FontSize',30);
    % Now add in second axis on the right:
    ax3 = axes('Position',get(ax2,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
    % Set units for the new axis using baseline values
    axis_limits = xlim(ax2);
    axis_limits_factor =[1 1]+axis_limits/100;
    new_limits =axis_limits_factor*baseline(i);
    xlim(ax3, new_limits);
    ylim(ax3, [0 1]);
    % Set width of axes
    % set(ax1, 'linewidth', 2);
    set(ax2, 'linewidth', 2);
    set(ax3, 'linewidth', 2);
    
    % make background transparent
    % set(ax1,'color','none');
    set(ax2,'color','none');
    set(ax3,'color','none');
    
    % remove new y axis
    set(ax3,'ycolor',get(fig,'color'),'ytick',[]);
    
    % remove box
    set(ax2,'Box','off');
    set(ax3,'Box','off');
    
    % label second axis
    actual_value_label = sprintf('%s %s, BCL = %0.f ms',FigTitleString, units, BCL(i));
    xlabel(ax3, actual_value_label);
    ylabel(ax2, 'Frequency')
    xlabh = get(ax3,'XLabel');
    set(xlabh,'Position',get(xlabh,'Position') + [0 -0.001 0]);
    
    % plot the 'extreme' point, ie 60 60 60 60 -60.
    hold on
    plot(ax3, results_60s(i), 0, 'kx', 'MarkerSize', 20, 'LineWidth', 4);
    hold off
    % set font size for all text
    hAll = findall(gcf);
    for idx = 1 : length(hAll)
        try
            set(hAll(idx),'fontsize',30);
        catch
            % nothing
        end
    end
             
    sname=sprintf('%s_hist_BCL%.0f',FileNameString,BCL(i));
    %print('-depsc',sname);
    export_fig([sname,'.png'], '-png', '-transparent')
    close all;
end
end

function PlotHistogramsActualValuesOnly(BCL, NonFailing, FailingWithVar, baseline, results_60s, NonFailingAltIndices, FailingWithVarAltIndices, FileNameString, FigTitleString, units)
% note: don't need explanatory data values for this!
exp_data=zeros(16385,11);

NonFailing = ((NonFailing/100) + 1);
FailingWithVar = ((FailingWithVar/100) + 1); 

% find min and max

min_vector = zeros(1,11);
max_vector = zeros(1,11);
max_hist_vector = zeros(1,11);
for i=1:length(BCL)
    [~,NonFailingNoAlt] = remove_alternans(exp_data, NonFailing(:,i), NonFailingAltIndices(:,1), NonFailingAltIndices(:,2), i);
    [~,FailingWithVarNoAlt] = remove_alternans(exp_data, FailingWithVar(:,i), FailingWithVarAltIndices(:,1), FailingWithVarAltIndices(:,2), i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I_NF = find(NonFailingNoAlt==-99); %(-10000/100) +1
    I_F = find(FailingWithVarNoAlt==-99); %(-10000/100) +1
    NonFailingNoAlt(I_NF) = [];
    FailingWithVarNoAlt(I_F) = [];
    min_vector(i) = min(min(NonFailingNoAlt*baseline(i)), min(FailingWithVarNoAlt*baseline(i)));
    max_vector(i) = max(max(NonFailingNoAlt*baseline(i)), max(FailingWithVarNoAlt*baseline(i)));

    % Determine bin width from smaller histogram
    range_F = max(FailingWithVarNoAlt*baseline(i)) - min(FailingWithVarNoAlt*baseline(i));
    range_NF = max(NonFailingNoAlt*baseline(i)) - min(NonFailingNoAlt*baseline(i));
    [range_for_binning,popn] = min([range_F range_NF]);
    if popn(1)==1 % considering case they're equal too...
        num_bins_failing = fix(0.3*sqrt(length(FailingWithVarNoAlt)));
        edges_failing = linspace(min(FailingWithVarNoAlt*baseline(i)), max(FailingWithVarNoAlt*baseline(i)), num_bins_failing);
        width_bins = edges_failing(2)-edges_failing(1);
        num_bins_nonfailing = fix(range_NF/width_bins);
        edges_nonfailing = linspace(min(NonFailingNoAlt*baseline(i)), max(NonFailingNoAlt*baseline(i)), num_bins_nonfailing);
    else
        num_bins_nonfailing = fix(0.3*sqrt(length(NonFailingNoAlt)));
        edges_nonfailing = linspace(min(NonFailingNoAlt*baseline(i)), max(NonFailingNoAlt*baseline(i)), num_bins_nonfailing);
        width_bins = edges_nonfailing(2)-edges_nonfailing(1);
        num_bins_failing = fix(range_F/width_bins);
        edges_failing = linspace(min(FailingWithVarNoAlt*baseline(i)), max(FailingWithVarNoAlt*baseline(i)), num_bins_failing);
    end
    [n1,x1]=histc(FailingWithVarNoAlt*baseline(i),edges_failing);
    [n2,x2]=hist(NonFailingNoAlt*baseline(i),edges_nonfailing);
    max_hist_vector(i) = max(max(n1),max(n2));
    
    % write statistics:
    statsFileName = ['Statistics_' FileNameString '_' num2str(BCL(i)) '.dat'];
    
    fp=fopen(statsFileName,'w');
    fprintf(fp, 'Non-Failing\n');
    fprintf(fp, 'mean\n');
    fprintf(fp, '%12.4e\n', mean(NonFailingNoAlt*baseline(i)));
    fprintf(fp, 'median\n');
    fprintf(fp, '%12.4e\n', median(NonFailingNoAlt*baseline(i)));
    fprintf(fp, 'standard deviation\n');
    fprintf(fp, '%12.4e\n', std(NonFailingNoAlt*baseline(i)));
    fprintf(fp, 'Inter-quartile range\n');
    fprintf(fp, '%12.4e\n', iqr(NonFailingNoAlt*baseline(i)));
    fprintf(fp, 'minimum\n');
    fprintf(fp, '%12.4e\n', min(NonFailingNoAlt*baseline(i)));
    fprintf(fp, 'maximum\n');
    fprintf(fp, '%12.4e\n', max(NonFailingNoAlt*baseline(i)));
    fprintf(fp, '\n');
    fprintf(fp, 'Failing\n');
    fprintf(fp, 'mean\n');
    fprintf(fp, '%12.4e\n', mean(FailingWithVarNoAlt*baseline(i)));
    fprintf(fp, 'median\n');
    fprintf(fp, '%12.4e\n', median(FailingWithVarNoAlt*baseline(i)));
    fprintf(fp, 'standard deviation\n');
    fprintf(fp, '%12.4e\n', std(FailingWithVarNoAlt*baseline(i)));
    fprintf(fp, 'Inter-quartile range\n');
    fprintf(fp, '%12.4e\n', iqr(FailingWithVarNoAlt*baseline(i)));
    fprintf(fp, 'minimum\n');
    fprintf(fp, '%12.4e\n', min(FailingWithVarNoAlt*baseline(i)));
    fprintf(fp, 'maximum\n');
    fprintf(fp, '%12.4e\n', max(FailingWithVarNoAlt*baseline(i)));
    fprintf(fp, '\n');
    perc_mean_change = 100*(mean(FailingWithVarNoAlt)-mean(NonFailingNoAlt))/mean(NonFailingNoAlt);
    perc_median_change = 100*(median(FailingWithVarNoAlt)-median(NonFailingNoAlt))/median(NonFailingNoAlt);
    fprintf(fp, '%% change in mean\n');
    fprintf(fp, '%12.4e\n', perc_mean_change);
    fprintf(fp, '%% change in median\n');
    fprintf(fp, '%12.4e\n', perc_median_change);

    fclose(fp);
  
    
end

data_min = min(min_vector);
data_max = max(max_vector);
hist_max = max(max_hist_vector);


for i=1:length(BCL)
    [~,NonFailingNoAlt] = remove_alternans(exp_data, NonFailing(:,i), NonFailingAltIndices(:,1), NonFailingAltIndices(:,2), i);
    [~,FailingWithVarNoAlt] = remove_alternans(exp_data, FailingWithVar(:,i), FailingWithVarAltIndices(:,1), FailingWithVarAltIndices(:,2), i);
    % now alternans is removed, we remove any points which had a '-1' reported
    % against them in the code.
    I_NF = find(NonFailingNoAlt==-99); %(-10000/100) +1
    I_F = find(FailingWithVarNoAlt==-99); %(-10000/100) +1
    NonFailingNoAlt(I_NF) = [];
    FailingWithVarNoAlt(I_F) = [];
    % nbins = fix(sqrt(length(NonFailingNoAlt)));
    
    % Determine bin width from smaller histogram
    range_F = max(FailingWithVarNoAlt*baseline(i)) - min(FailingWithVarNoAlt*baseline(i));
    range_NF = max(NonFailingNoAlt*baseline(i)) - min(NonFailingNoAlt*baseline(i));
    
    [range_for_binning,popn] = min([range_F range_NF]);
    if popn(1)==1 % considering case they're equal too...
        num_bins_failing = fix(0.3*sqrt(length(FailingWithVarNoAlt)));
        edges_failing = linspace(min(FailingWithVarNoAlt*baseline(i)), max(FailingWithVarNoAlt*baseline(i)), num_bins_failing);
        width_bins = edges_failing(2)-edges_failing(1);
        num_bins_nonfailing = fix(range_NF/width_bins);
        edges_nonfailing = linspace(min(NonFailingNoAlt*baseline(i)), max(NonFailingNoAlt*baseline(i)), num_bins_nonfailing);
    else
        num_bins_nonfailing = fix(0.3*sqrt(length(NonFailingNoAlt)));
        edges_nonfailing = linspace(min(NonFailingNoAlt*baseline(i)), max(NonFailingNoAlt*baseline(i)), num_bins_nonfailing);
        width_bins = edges_nonfailing(2)-edges_nonfailing(1);
        num_bins_failing = fix(range_F/width_bins);
        edges_failing = linspace(min(FailingWithVarNoAlt*baseline(i)), max(FailingWithVarNoAlt*baseline(i)), num_bins_failing);
    end
    
    fig=figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]);
    hold on
    [n1,x1]=histc(FailingWithVarNoAlt*baseline(i),edges_failing);
    h1 =bar(edges_failing,n1,'histc');
    %set(h1,'FaceColor',0.55*[1 1 1],'EdgeColor','none','facealpha',1.0, 'edgealpha', 0.5);
    h = findobj(gca,'Type','patch');
    %set(h,'FaceColor','r','EdgeColor','r','facealpha',0.75, 'edgealpha', 0.75);
    set(h,'FaceColor','r','EdgeColor','none','facealpha',0.75, 'edgealpha', 0.75);
    [n2,x2]=hist(NonFailingNoAlt*baseline(i),edges_nonfailing);
    h2 =bar(edges_nonfailing,n2,'histc');
    %set(h2,'facecolor',0*[1 1 1],'facealpha',0.75, 'EdgeColor', 'none', 'edgealpha', 0.5);
    h = findobj(gca,'Type','patch');
    %set(h,'facealpha',0.75, 'EdgeColor', 'b', 'edgealpha', 0.75);
    set(h,'facealpha',0.75, 'EdgeColor', 'none', 'edgealpha', 0.75);
    hold off
    ax2 = gca;
    
    %current_position=get(ax2,'Position');
    %set(ax2, 'Position',current_position+[0 -0.015 0 -0.1]);
    
    %histleg=legend('HF','Normal');
    %set(histleg, 'Position', get(histleg, 'Position')+[-0.025 -0.1 0 0]);
    %xlabel('% Change from baseline','FontSize',30);
     % Set width of axes
    % set(ax1, 'linewidth', 2);
    set(ax2, 'linewidth', 2);

    % make background transparent
    % set(ax1,'color','none');
    set(ax2,'color','none');
  
    % remove box
    set(ax2,'Box','off');
    
    % set axes limits
    set(ax2, 'xlim',[data_min data_max], 'ylim', [0 hist_max]);
    
    %set(xlabel_handle, 'Position', get(xlabel_handle, 'Position')+[0 -12 0]);
    %ylabel(ax2, 'Frequency')
    
    % set font size for all text
    hAll = findall(gcf);
    for idx = 1 : length(hAll)
        try
            set(hAll(idx),'fontsize',26);
        catch
            % nothing
        end
    end
    % label second axis
    %actual_value_label = sprintf('%s %s',FigTitleString, units);
    %xlabel_handle=xlabel(ax2, actual_value_label, 'FontSize', 30);
    
    ax3=axes('Position', get(ax2, 'Position'));
    set(ax3,'xlim', get(ax2, 'xlim'));
    set(ax3,'ylim', get(ax2, 'ylim'));
    set(ax3, 'color', 'none');
    set(ax3, 'xtick', []);
    set(ax3, 'ytick', []);
    
    hold on
    % plot the 'extreme' point, ie 60 60 60 60 -60.
    if results_60s(i)~=-1
        plot(ax3,results_60s(i), 0,  'kx', 'MarkerSize', 20, 'LineWidth', 4);
    end
    % plot the 'baseline' point, ie 0 0 0 0 0
    plot(ax3,baseline(i), 0, 'ko', 'MarkerSize', 20, 'LineWidth', 4);
    hold off
             
    sname=sprintf('%s_hist_BCL%.0f',FileNameString,BCL(i));
    %print('-depsc',sname);
    %export_fig([sname,'.eps'], '-eps', '-transparent')
    pause(10)
    export_fig([sname,'.png'], '-png', '-transparent')
    pause(10);
    close all;
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
