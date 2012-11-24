function exp_des_uniform
%
% Mask for parameters to include in the sensitivity 
%
% Parameters vector and mask
%
% ICaL : 0
% Tauf : 1
% Tauf2: 2
% GKr  : 3
% GKs  : 4
% TauXs: 5
% GK1  : 6
% GNaK : 7
% GNaCa: 8
% Gto  : 9
% Iup  : 10
%
% Setting the mask for the modifiable parameters
%
% This mask will change
% ICaL, GKr, GKs, GK1, GNaCa, Gto, Iup
%
mask=[0 3 8 9 10];
%
% parameters that are down-, up, and non-rgulated in the experiment
% indexing correspond to the mask vector
%
downreg = [1 2 4 5];
upreg = [3];
nonreg = [];
% downreg = [];
% upreg = [];
% nonreg = 1:7;
%
% percentage for down-, up-, and none- regulation.
%
perc=[-0.6 0.6 0.3];
% perc=[0.0 0.0 0.3];
blck = 656; % make divisible by 8 to run on HAL.
tag = '3060nonvar'; % for experiment file name Id
gen_exp_des_uniform(mask,downreg,upreg,nonreg,perc,blck,tag)
return

function gen_exp_des_uniform(mask,downreg,upreg,nonreg,perc,blck,tag)
%
% writes a 3 level full factorial experimental design file
% for the simulations according to the size of mask
%
% perc is the percentage of variation, i.e. 0.6 -> 60%
% blck in how many blocks we want to split the experiment design for
% computation
%
nvar = length(mask);
%
% generate nlevel^7 random parameter sets in the range [0,1];
%
nlevel=4;
unif_dist = rand(nlevel^7,nvar);
%
% setting the levels for the experiment accoridng to the perc vector
%
% downregulation
%
if (length(downreg)>0)
  unif_dist(:,downreg) = 1+perc(1)*unif_dist(:,downreg);
end
%
% upregulation
%
if (length(upreg)>0)
  unif_dist(:,upreg) = 1+perc(2)*unif_dist(:,upreg);
end
%
% noneregulation (central design)
%
if (length(nonreg)>0)
  unif_dist(:,nonreg) = 1+perc(3)*(2*unif_dist(:,nonreg)-1.0);
end
%
n=length(unif_dist(:,1));
cont_exp=1;
fname=sprintf('exp_design_uniform_%s_%dvar.dat',tag,nvar);
fp=fopen(fname,'w');
fprintf(fp,'%6d\n',n);
fprintf(fp,'%6d ',nvar);
fprintf(fp,'%6d ',mask);
fprintf(fp,'\n');
for j=1:n
    fprintf(fp,'%6.0f ',1.0*cont_exp);
%    for k=nvar:-1:1
%        fprintf(fp,'%6.4f ',unif_dist(cont_exp,k));
%    end
    fprintf(fp,'%6.4f ',unif_dist(cont_exp,:));
    fprintf(fp,'\n');
    cont_exp = cont_exp + 1;
end
if blck>0
    Nmin = fix(n/blck);
    Npb = Nmin*ones(blck,1);
    nxtra = mod(n,blck);
    Npb(1:nxtra)=Npb(1:nxtra)+1;
    cont_exp=1;
    jend=0;
    for i=1:blck
        jini = jend+1;
        jend = jend+Npb(i);
        fname=sprintf('exp_design_uniform_%s_%d.dat',tag,i);
        fp=fopen(fname,'w');
        fprintf(fp,'%6d\n',jend-jini+1);
        fprintf(fp,'%6d ',nvar);
        fprintf(fp,'%6d ',mask);
        fprintf(fp,'\n');
        for j=jini:jend
            fprintf(fp,'%6.0f ',1.0*cont_exp);
%            for k=nvar:-1:1
%                fprintf(fp,'%6.4f ',unif_dist(cont_exp,k));
%            end
            fprintf(fp,'%6.4f ',unif_dist(cont_exp,:));
            fprintf(fp,'\n');
            cont_exp = cont_exp + 1;
        end
        if(jend>n)
            jend=n;
        end
    end
end
return
