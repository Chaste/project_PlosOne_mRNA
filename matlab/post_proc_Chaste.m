function post_proc
%
% Computes the following biomarkers:
%  i) APD80
% ii) APD30/APD80 ratio
%iii) CaiTD
% iv) CaiTD30/CaiTD80 ratio
%  v) AP Cai delay
% vi) Systolic Cai
%
% Data is normalized with respect to the baseline value
%
% postprocessing program
%
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
% 'dirname': is the directory where the simulated data is located. 
% 'fname'  : contains the experiment design used in the simulations. 
%            fname has to be located in the same directory
%
% dirname='../results/exp_des_uniform_4_30.HAL.260312/';
% fname='exp_design_uniform_30_7var.dat';
% dirname='../results/exp_des_uniform_4_3060.HAL.280312/';
% fname='exp_design_uniform_3060_7var.dat';
% dirname='../results/exp_des_uniform_4_3060nonvar.HAL.290312/';
% fname='exp_design_uniform_3060nonvar_5var.dat';
% dirname='../results/exp_design_30_surf/';
% fname='exp_design_30_surf_2var.dat';
% dirname='../results/exp_design_60_surf/';
% fname='exp_design_60_surf_2var.dat';
%dirname='../results/OHara2011_endo_uniform_30/';
%fname='exp_design_uniform_30_7var.dat';
%dirname='../results/OHara2011_endo_uniform_3060/';
%fname='exp_design_uniform_3060_7var.dat';
%dirname='../results/OHara2011_endo_uniform_3060nonvar/';
%fname='exp_design_uniform_3060nonvar_5var.dat';
%dirname='../results/OHara2011_endo_30_surf/';
%fname='exp_design_30_surf_2var.dat';
%dirname='../results/OHara2011_endo_60_surf/';
%fname='exp_design_60_surf_2var.dat';
dirname='../results/OHara2011_endo_extremecase/';
fname='exp_design_extremecase.dat';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% we start computations
% 
cd(dirname);
[sucess,message,messageid]=mkdir('post_data');
if(sucess==0)
    fprintf('dir post_data. %s. Matlab error: %s\n',message,messageid);
    return;
elseif(isempty(message))
    fprintf('dir post_data. %s. Matlab error: %s\n',message,messageid);
end
%fac = 0.5/per;
fp = fopen(fname,'r');
Nexp = fscanf(fp,'%d',1);      % Total number of sample points
Nvar = fscanf(fp,'%d',1);      % Number of variables
mask = fscanf(fp,'%d',Nvar);   % Mask with model parameters
ff = zeros(Nexp,Nvar);        % reading the file with the
for i=1:Nexp                   % experimental design
    num = fscanf(fp,'%d',1);
    for j=1:Nvar
        ff(i,j)=fscanf(fp,'%f',1);
    end
end
%ff = (ff-1.0);
fclose(fp);
nexp_base = find(prod(ff,2)==1);
%
% directory where post-process biomarkers will be writen
%
rsa_dir=sprintf('%s/%spost_data',curr_dir,dirname);
%
write_info(Nexp,Nvar,mask,fname,dirname);
%write_table_parameters(mask,curr_dir);
%
% COMPUTING BIOMARKERS
%
%
% Dynamic protocol
%
fprintf('postprocessing dynamic protocol\n');


strArray = java_array('java.lang.String', 11);
strArray(1) = java.lang.String(' APD80 ');
strArray(2) = java.lang.String('APD3080');
strArray(3) = java.lang.String('CaiTD80');
strArray(4) = java.lang.String('Cai3080');
strArray(5) = java.lang.String('DAPCaiT');
strArray(6) = java.lang.String('SystCai');
strArray(7) = java.lang.String(' APAlt ');
strArray(8) = java.lang.String(' Ca80Alt ');
strArray(9) = java.lang.String('2Cai3080 ');
strArray(10) = java.lang.String('CaiTD30 ');
strArray(11) = java.lang.String('Ca30Alt ');
tagSS = cell(strArray);

[BCL,APD80,APD3080,CaiTD80,Cai3080,DAPDCaiD,SysCai,AP_Alternans, CaiTD80_Alternans, Cai3080_2, CaiTD30, CaiTD30_Alternans]=post_protocol(Nexp);

% AP_Alternans

%nexp_base = ceil(Nexp/2);
basevalues = [APD80(nexp_base,:); APD3080(nexp_base,:); ...
              CaiTD80(nexp_base,:); Cai3080(nexp_base,:); ...
              DAPDCaiD(nexp_base,:);SysCai(nexp_base,:); ...
              AP_Alternans(nexp_base,:);CaiTD80_Alternans(nexp_base,:);...
              Cai3080_2(nexp_base,:); CaiTD30(nexp_base,:);...
              CaiTD30_Alternans(nexp_base,:)];
%
% normalizing the output with respect to baseline
%
nBCL = length(BCL);
for i=1:nBCL
    for n=1:Nexp
        if APD80(n,i)~=-1
            APD80(n,i) = 100*(APD80(n,i)-basevalues(1,i))/basevalues(1,i);
            APD3080(n,i) = 100*(APD3080(n,i)-basevalues(2,i))/basevalues(2,i);
        else 
            APD80(n,i)=-10000;
            APD3080(n,i)=-10000;
        end
        if CaiTD80(n,i)~=-1
            CaiTD80(n,i) = 100*(CaiTD80(n,i)-basevalues(3,i))/basevalues(3,i);
            Cai3080(n,i) = 100*(Cai3080(n,i)-basevalues(4,i))/basevalues(4,i);
            Cai3080_2(n,i) = 100*(Cai3080_2(n,i)-basevalues(9,i))/basevalues(9,i);
            CaiTD30(n,i) = 100*(CaiTD30(n,i)-basevalues(10,i))/basevalues(10,i);
        else
            CaiTD80(n,i) = -10000;
            Cai3080(n,i) = -10000;
            Cai3080_2(n,i) = -10000;
            CaiTD30(n,i) = -10000;
        end
        if DAPDCaiD(n,i)~=-1
            DAPDCaiD(n,i) = 100*(DAPDCaiD(n,i)-basevalues(5,i))/basevalues(5,i);
        else
            DAPDCaiD(n,i) = -10000;
        end
        if SysCai(n,i)~=-1
            SysCai(n,i) = 100*(SysCai(n,i)-basevalues(6,i))/basevalues(6,i);
        else
            SysCai(n,i) = -10000;
        end
    end
end

%
% writing baseline values
%
write_value_Baseline_Biomarkers(BCL,basevalues,rsa_dir);
%
%
% writing scaled posprocessed values
%
write_scaled_bio(rsa_dir,'APD80',tagSS{1},Nexp,Nexp,APD80);
write_scaled_bio(rsa_dir,'APD3080',tagSS{2},Nexp,Nexp,APD3080);
write_scaled_bio(rsa_dir,'CaiTD80',tagSS{3},Nexp,Nexp,CaiTD80);
write_scaled_bio(rsa_dir,'Cai3080',tagSS{4},Nexp,Nexp,Cai3080);
write_scaled_bio(rsa_dir,'DAPCaiT',tagSS{5},Nexp,Nexp,DAPDCaiD);
write_scaled_bio(rsa_dir,'SysCai',tagSS{6},Nexp,Nexp,SysCai);
write_scaled_bio(rsa_dir,'APAlt',tagSS{7},Nexp,Nexp,AP_Alternans);
write_scaled_bio(rsa_dir,'CaAlt',tagSS{8},Nexp,Nexp,CaiTD80_Alternans);
write_scaled_bio(rsa_dir,'2Cai3080',tagSS{9},Nexp,Nexp,Cai3080_2);
write_scaled_bio(rsa_dir,'CaiTD30',tagSS{10},Nexp,Nexp,CaiTD30);
write_scaled_bio(rsa_dir,'Ca30Alt',tagSS{11},Nexp,Nexp,CaiTD30_Alternans);
%
% writing the experiment design
%
write_exp_design(rsa_dir,mask,ff)


% Find the files which did not record biomarkers
[apd30nr, apd80nr, caitd30nr, caitd80nr, maxcaitnr, delaynr]=find_not_recorded(Nexp, BCL);

% Write biomarker files not recorded to disk
write_not_recorded_list(rsa_dir, 'APD30', apd30nr,BCL);
write_not_recorded_list(rsa_dir, 'APD80', apd80nr,BCL);
write_not_recorded_list(rsa_dir, 'CaiTD30', caitd30nr,BCL);
write_not_recorded_list(rsa_dir, 'CaiTD80', caitd80nr,BCL);
write_not_recorded_list(rsa_dir, 'delay', delaynr,BCL);
write_not_recorded_list(rsa_dir, 'maxcait', maxcaitnr,BCL);

cd(curr_dir);
end
%
%  Additional functions
%
function write_info(Nexp,Nvar,mask,fname,dirname)
%
%
fprintf('Processing experiment: %s\n',fname);
fprintf('Source directory     : %s\n',dirname);
fprintf('Number of experiments: %d\n',Nexp);
fprintf('Number of variables  : %d\n',Nvar);
fprintf('Variables in the experiment:\n');
for i=1:Nvar
    switch mask(i)
        case 0
            fprintf('GCaL ');
        case 1
            fprintf('tau_f ');
        case 2
            fprintf('tau_f2 ');
        case 3
            fprintf('GKr ');
        case 4
            fprintf('GKs ');
        case 5
            fprintf('tau_Xs ');
        case 6
            fprintf('GK1 ');
        case 7
            fprintf('GNaK ');
        case 8
            fprintf('GNCx ');
        case 9
            fprintf('Gto ');
        case 10
            fprintf('Gup ');
    end
end
fprintf('\n\n');
end

%
function write_value_Baseline_Biomarkers(BCL,baseline,rsa_dir)
%
%
strArray = java_array('java.lang.String', 12);
strArray(1) = java.lang.String('     BCL    ');
strArray(2) = java.lang.String('    APD80   ');
strArray(3) = java.lang.String('   APD3080  ');
strArray(4) = java.lang.String('    CaiTD   ');
strArray(5) = java.lang.String('   Cai3080  ');
strArray(6) = java.lang.String('   DAPCaiD  ');
strArray(7) = java.lang.String('   SysCai   ');
strArray(8) = java.lang.String(' AP Alternans');
strArray(9) = java.lang.String(' Ca Alternans');
strArray(10) = java.lang.String(' Cai3080 pb');
strArray(11) = java.lang.String(' CaiTD30    ');
strArray(12) = java.lang.String(' Cai30Alt   ');
tagDYN = cell(strArray);
filename = sprintf('%s/baselineBio.dat',rsa_dir);
fp=fopen(filename,'w');
fprintf(fp,'%%');
for i=1:length(tagDYN)
 fprintf(fp,'%s ',tagDYN{i});
end
fprintf(fp,'\n');
for i=1:length(BCL)
    fprintf(fp,'%12.4e ',BCL(i));
    fprintf(fp,'%12.4e ',baseline(:,i));
    fprintf(fp,'\n');
end
fclose(fp);
end
%
function write_exp_design(rsa_dir,mask,ff)
%
filename = sprintf('%s/exp_design.dat',rsa_dir);
fp=fopen(filename,'w');
Nvar = length(mask);
Nexp = length(ff(:,1));
fprintf(fp,'%%');
for i=1:Nvar
    switch mask(i)
        case 0
            fprintf(fp,' GCaL  ');
        case 1
            fprintf(fp,'  tau_f ');
        case 2
            fprintf(fp,' tau_f2 ');
        case 3
            fprintf(fp,'   GKr  ');
        case 4
            fprintf(fp,'   GKs  ');
        case 5
            fprintf(fp,' tau_Xs ');
        case 6
            fprintf(fp,'   GK1  ');
        case 7
            fprintf(fp,'  GNaK  ');
        case 8
            fprintf(fp,'  GNCx  ');
        case 9
            fprintf(fp,'   Gto  ');
        case 10
            fprintf(fp,'   Gup  ');
    end
end
fprintf(fp,'\n');
fprintf(fp,'%d ',mask);
fprintf(fp,'\n');
for i=1:Nexp
    fprintf(fp,'%7.4f ',ff(i,:));
    fprintf(fp,'\n');
end
fclose(fp);
end


function write_scaled_bio(rsa_dir,fname,tags,n_exp_fac,Nexp,Vec)
%
% writing scaled posprocessed values
%
filename = sprintf('%s/%s.dat',rsa_dir,fname);
fp=fopen(filename,'w');
fprintf(fp,'%%');
for i=1:length(tags(:,1))
 fprintf(fp,'%s ',tags(i,:));
end
fprintf(fp,'\n');
for i=1:Nexp
    fprintf(fp,'%12.4e ',Vec(i,:)/100);
    fprintf(fp,'\n');
end
fclose(fp);    
end

function [apd30nr, apd80nr, caitd30nr, caitd80nr, maxcaitnr, delaynr]=find_not_recorded(Nexp, BCL)

for j=1:length(BCL)
    apd30nr{j} = [];
    apd80nr{j} = [];
    caitd30nr{j} = [];
    caitd80nr{j} = [];
    maxcaitnr{j} = [];
    delaynr{j} = [];
end

for i=1:Nexp
    fname = sprintf('data/DYN_%d.out', i);
    a=load(fname);
    for j=1:length(BCL)
        apd30list_current = apd30nr{j};
        apd80list_current = apd80nr{j};
        caitd30list_current = caitd30nr{j};
        caitd80list_current = caitd80nr{j};
        delaylist_current = delaynr{j};
        maxcaitlist_current = maxcaitnr{j};

        if (a(j,2) == -1 || a(j,3) == -1)
            apd30list_current = [apd30list_current; i];
            apd30nr{j} = apd30list_current;
        end
        if (a(j,4) == -1 || a(j,5) == -1)
            apd80list_current = [apd80list_current; i];
            apd80nr{j} = apd80list_current;
        end
        if (a(j,6) == -1 || a(j,7) == -1)
            caitd30list_current = [caitd30list_current; i];
            caitd30nr{j} = caitd30list_current;
        end
        if (a(j,8) == -1 || a(j,9) == -1)
            caitd80list_current = [caitd80list_current; i];
            caitd80nr{j} = caitd80list_current;
        end
        if (a(j,10) == -1 || a(j,11) == -1)
            delaylist_current = [delaylist_current; i];
            delaynr{j} = delaylist_current;
        end
        if (a(j,12) == -1 || a(j,13) == -1)
            maxcaitlist_current = [maxcaitlist_current; i];
            maxcaitnr{j} = maxcaitlist_current;
        end
    end
    
end
 
end


function write_not_recorded_list(rsa_dir, biomarker_string, biomarker_cell,BCL)

for i=1:length(BCL)
    exp_list = biomarker_cell{i};
    if ~isempty(exp_list)
        filename = sprintf('%s/%sNotRecorded%d.dat',rsa_dir,biomarker_string, BCL(i));
        listfile = fopen(filename,'w');
        for j=1:length(exp_list)
            fprintf(listfile, '%d\n', exp_list(j));
        end
    end
end    

end

%
% Biomarker Postprocessing functions
%
function [BCL,APD80,APD3080,CaiTD80,CaiD3080_1,DAPDCaiD,SysCai, AP_Alternans, CaiTD80_Alternans, CaiD3080_2, CaiTD30, CaiTD30_Alternans]=post_protocol(Nexp)
%
% postprocessing dynamic protocol
%
alternans_threshold=5;

fname = sprintf('data/DYN_1.out');
a=load(fname);
nBCL = length(a(:,1));
BCL = a(:,1)';
APD80 = zeros(Nexp,nBCL);
APD3080 = zeros(Nexp,nBCL);
CaiTD80 = zeros(Nexp,nBCL);
CaiD3080_1 = zeros(Nexp,nBCL);
CaiD3080_2 = zeros(Nexp,nBCL);
DAPDCaiD = zeros(Nexp,nBCL);
SysCai = zeros(Nexp,nBCL);
AP_Alternans = zeros(Nexp,nBCL);
CaiTD30 = zeros(Nexp,nBCL);
CaiTD30_Alternans = zeros(Nexp,nBCL);
APD80(1,:)    = a(:,5)';
APD3080(1,:)  = (a(:,3)./a(:,5))';
CaiTD80(1,:)     = a(:,9)';
CaiD3080_1(1,:) = (a(:,7)./a(:,9))';
CaiD3080_2(1,:) = (a(:,6)./a(:,8))';
DAPDCaiD(1,:) = a(:,11)';
SysCai(1,:) = a(:,13)';
APD80_difference = a(:,5)'-a(:,4)';
AlternansLong=(APD80_difference>alternans_threshold);
AlternansShort=(APD80_difference<-alternans_threshold);
AP_Alternans(1,:)=AlternansLong+2*AlternansShort;
CaiTD80_difference = a(:,9)'-a(:,8)';
CaiTD80AlternansLong=(CaiTD80_difference>alternans_threshold);
CaiTD80AlternansShort=(CaiTD80_difference<-alternans_threshold);
CaiTD80_Alternans(1,:)=CaiTD80AlternansLong+2*CaiTD80AlternansShort;
% CaiTD30
CaiTD30(1,:) = a(:,7)';
CaiTD30_difference = a(:,7)'-a(:,6)';
CaiTD30AlternansLong=(CaiTD30_difference>alternans_threshold);
CaiTD30AlternansShort=(CaiTD30_difference<-alternans_threshold);
CaiTD30_Alternans(1,:)=CaiTD30AlternansLong+2*CaiTD30AlternansShort;

for j=1:length(BCL)
    if a(j,4) == -1 || a(j,5) == -1
        AP_Alternans(1,j) = 3;
    end
    if a(j,6) == -1 || a(j,7) == -1
        CaiTD30_Alternans(1,j) = 3;
    end
    if a(j,8) == -1 || a(j,9) == -1
        CaiTD80_Alternans(1,j) = 3;
    end
end

for i=2:Nexp
    fname = sprintf('data/DYN_%d.out',i);
    a=load(fname);
    APD80(i,:)    = a(:,5)';
    APD3080(i,:)  = (a(:,3)./a(:,5))';
    CaiTD80(i,:)     = a(:,9)';
    CaiD3080_1(i,:) = (a(:,7)./a(:,9))';
    CaiD3080_2(i,:) = (a(:,6)./a(:,8))';
    DAPDCaiD(i,:) = a(:,11)';
    SysCai(i,:) = a(:,13)';
    APD80_difference = a(:,5)'-a(:,4)';
    AlternansLong=(APD80_difference>alternans_threshold);
    AlternansShort=(APD80_difference<-alternans_threshold);
    AP_Alternans(i,:)=AlternansLong+2*AlternansShort;
    CaiTD80_difference = a(:,9)'-a(:,8)';
    CaiTD80AlternansLong=(CaiTD80_difference>alternans_threshold);
    CaiTD80AlternansShort=(CaiTD80_difference<-alternans_threshold);
    CaiTD80_Alternans(i,:)=CaiTD80AlternansLong+2*CaiTD80AlternansShort;
    %CaiT30
    CaiTD30(i,:)= a(:,7)';
    CaiTD30_difference = a(:,7)'-a(:,6)';
    CaiTD30AlternansLong=(CaiTD30_difference>alternans_threshold);
    CaiTD30AlternansShort=(CaiTD30_difference<-alternans_threshold);  
    CaiTD30_Alternans(i,:)=CaiTD30AlternansLong+2*CaiTD30AlternansShort;
    for j=1:length(BCL)
    if a(j,4) == -1 || a(j,5) == -1
        AP_Alternans(i,j) = 3;
    end
    if a(j,6) == -1 || a(j,7) == -1
        CaiTD30_Alternans(i,j) = 3;
    end
    if a(j,8) == -1 || a(j,9) == -1
        CaiTD80_Alternans(i,j) = 3;
    end
end
end
end

