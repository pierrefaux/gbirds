% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Pierre FAUX (pierrefaux@gmail.com), 2019
% GBIRDS Version 1.2, last update on 2019-10-04


% This script uses the outputs of GBIRDS, namely
%   - sel3.012
%   - sel3.012.markers
%   - sel3.012.indivs
%
% The script requires the following functions:
%   - GgbsGBirdS.m (Adapted from Dodds et al. 2015, BMC Genomics)
%   - Fst3GBirdS.m
%   - getFstGBirdS.m
%   - plotKGBirdS.m
%
% USAGE:
%   Run this function in a Matlab/Octave session in the directory with
%   GBirdS outputs (sel3.012*)
%

function FinalSteps

% 1. DATA LOADING
gc=load('sel3.012');    % genotype calls of sel III
mark=dlmread('sel3.012.markers',' ',1,19);  % marker information
    cr=mark(:,1);   % call rate
[ids,pop]=textread('sel3.012.indiv','%s%d');
nind=length(ids);
ad=load('sel3.012.adepth');
v=(1:2:nind*2-1);
da1=ad(v,:);
da2=ad(v+1,:);

% 2. COMPUTATION OF Fst
gc1=gc(:,cr==1); % consider only the markers with 100% calls
getFstGBirdS(gc1,pop);

% 3. COMPUTATION OF KGD MATRIX
G=GgbsGBirdS(da1,da2,2,Inf);

% 4. PCA
[S,L]=pcacov(G);

% 5. PLOT OF PC2 vs. PC1
plotKGBirdS(S,L,ids,'',3);


% FUNCTIONS
function getFstGBirdS(GC2, pops)
[n,m]=size(GC2);
np=numel(unique(pops));
fprintf ('%d indv x %d markers\n', n, m);
for p1=1:np-1
	for p2=p1+1:np
		F=zeros(m,1);
		fid=fopen(strcat('Fst3_pop',num2str(p1),'_vs_pop',num2str(p2)),'w');
		nfst=0;
		for i=1:m
			indsp1=find(GC2(:,i)>-1 & pops==p1);
			indsp2=find(GC2(:,i)>-1 & pops==p2);
			if numel(indsp1)>0 && numel(indsp2)>0
				inds=sort([indsp1;indsp2]);
				nm=length(inds);
				gc=GC2(inds,i);
				p=pops(inds);
				F(i)=Fst3GBirdS(gc,p);
				nfst=nfst+1;
			else
				nm=0;
				F(i)=NaN;
			end
			fprintf(fid,'%d %9.5f\n', nm, F(i));
		end
		mf=mean(F(isnan(F)==0));
		sf=std(F(isnan(F)==0));
		fprintf('POP%d vs. POP%d :: #Fst = %d avg Fst = %9.4f +- %9.4f\n', p1, p2, nfst, mf, sf);
		fclose(fid);
	end
end


function fst=Fst3GBirdS(gc,spop)
m=length(gc(1,:));
fst=zeros(m,1);
spops=unique(spop(spop>0));
r=numel(spops);
nis=zeros(r,1);
pis=zeros(r,1);
his=zeros(r,1);
for i=1:r
    nis(i)=numel(gc(spop==spops(i),1));
end
nim=sum(nis)/r;
nc=(r*nim-nis'*nis/(r*nim))/(r-1);
for j=1:m
    for i=1:r
        pis(i)=sum(gc(spop==spops(i),j))/(2*nis(i));
        his(i)=numel(find(gc(spop==spops(i),j)==1))/nis(i);
    end
    pm=nis'*pis/(r*nim);
    hm=nis'*his/(r*nim);
    s2=nis'*((pis-pm).^2)/((r-1)*nim);
    a=nim/nc*(s2-(pm*(1-pm)-(r-1)/r*s2-hm/4)/(nim-1));
    b=nim*(pm*(1-pm)-(r-1)/r*s2-(2*nim-1)*hm/(4*nim))/(nim-1);
    c=.5*hm;
    fst(j)=a/(a+b+c);
end


function plotKGBirdS(K,L,id,titre,nc)
n=length(id);
fprintf('#ids= %d\n', n);
u=char(id);
lab=unique(u(:,1:nc),'rows');
nl=length(lab(:,1));
fprintf('#labels= %d\n', nl);
labs=zeros(n,3);
for i=1:n
    u=char(id{i});
    for j=1:nl
        if strcmp(u(1:nc),lab(j,:))
            labs(i,:)=mycolors(j+3);
        end
    end
end
figure('Color','w','Position',[1 1 1440 900]);
title(titre,'FontSize',24,'FontWeight','Bold');
z1=0.005*(max(K(:,1))-min(K(:,1)));
z2=0.005*(max(K(:,2))-min(K(:,2)));
for i=1:n
    hold on;
    plot(K(i,1),K(i,2),'x','LineWidth',3,'Color',labs(i,:));
    text(K(i,1)+z1,K(i,2)+z2,id{i},'FontSize',14,'FontWeight','Bold', ...
        'Color',labs(i,:));
end
l1=100*L(1)/sum(L);
l2=100*L(2)/sum(L);
xlabel(strcat('1^{st} principal component: ',num2str(l1,'%9.2f%%')), ...
    'FontSize',16,'FontWeight','Bold');
ylabel(strcat('2^{nd} principal component: ',num2str(l2,'%9.2f%%')), ...
    'FontSize',16,'FontWeight','Bold');
set(gca,'FontSize',16,'FontWeight','Bold');
set(gca,'Box','Off');


function G=GgbsGBirdS(DA1,DA2,depthmin,depthmax)
gc=zeros(size(DA1))*NaN;
gc(DA1>0 & DA2>0)=1;    %A1/A2
gc(DA1>0 & DA2==0)=2;   %A1/A1
gc(DA1==0 & DA2>0)=0;   %A2/A2
[nind,nsnps]=size(gc);
depth=DA1+DA2;
freq=zeros(1,nsnps);
for i=1:nsnps
    a=sum(DA1(:,i));
    b=sum(DA2(:,i));
    freq(i)=a/(a+b);
end
usegeno=zeros(size(gc));
usegeno(isnan(gc)==0)=1;	% usegeno is an integer (not logical) mask
if (depthmin>1 || depthmax < Inf)
	gc(depth<depthmin)=NaN;
	gc(depth>depthmax)=NaN;
	depth(depth<depthmin)=0;
	depth(depth>depthmax)=0;
	usegeno(depth<depthmin)=0;
	usegeno(depth>depthmax)=0;
end
for i=1:nind
	gc(i,:)=gc(i,:)-2*freq;
end
gc(isnan(gc))=0;
P0=kron(ones(nind,1),freq);
P1=1-P0;
P0(usegeno==0)=0;
P1(usegeno==0)=0;
div0=2*P0*P1';
G=gc*gc'./div0;
gc(depth<2)=0;
P0(depth<2)=0;
P1(depth<2)=0;
depth2K=2.^-(depth);
for i=1:nind
	G(i,i)=1+sum(...
		(gc(i,:).^2-2*P0(i,:).*P1(i,:).* ...
		(1+2*depth2K(i,:)))./ ...
		(1-2*depth2K(i,:)))/div0(i,i);
end
fprintf('Mean self-relatedness: %9.6d\n', mean(diag(G)));
