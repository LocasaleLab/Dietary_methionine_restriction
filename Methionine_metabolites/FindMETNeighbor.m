%-------------------------------------------------------
% Define methionine-related and methionine-unrelated
% metabolites based on the human genome-scale metabolic 
% model Recon2
%-------------------------------------------------------
load Recon2.mat                                         % Load the metabolic network model

MET={'L-methionine'};                                   %Name of methionine in the model
[~,b]=ismember(model.metNames,MET);                     
METPos=find(b>0);                                       %Find all species named methionine (in different compartments)
Reliable={'4'};                                         
[~,b1]=ismember(model.rxnConfidenceScores,Reliable);    %Only consider reactions with highest confidence scores
Null={''};
[a,b2]=ismember(model.rxnECNumbers,Null);               %Only consider enzyme-catalyzed reactions
ECPos=find(b1>0 & b2==0);

SSub=model.S(:,ECPos);                                  %Remove low-confidence and non-enzymatic reactions
NeighborPos{1}=METPos;                                  %Start from the metabolite set with distance=0
HubPos=find(sum(abs(model.S'))>100);                    %Identify hub metabolites (cofactors, water, etc)
for n=1:6
    CensorPos=NeighborPos{n};                           %Work on the nth set to define the n+1th set
    NeighborPos{n+1}=[];
    for i=1:length(CensorPos)
        pos=findOtherEnd(CensorPos(i),SSub);
        pos=setdiff(pos,HubPos);                        %Remove hub metabolites
        pos=pos(:);
        NeighborPos{n+1}=union(NeighborPos{n+1},pos);
    end
    NeighborPos{n+1}=union(NeighborPos{n},NeighborPos{n+1});
end

%-------------------------------------------------------------------
% Map KEGG IDs to metabolite names used in the metabolomics dataset
%-------------------------------------------------------------------
K2M_KEGG=textread('K2M_KEGG.txt','%s','delimiter','\n');
K2M_Met=textread('K2M_Met.txt','%s','delimiter','\n');
KEGG_D3=unique(model.metKeggID(NeighborPos{4}));        %KEGG IDs of metabolites within 3 steps from methionine
KEGG_D4=unique(model.metKeggID(NeighborPos{5}));        %KEGG IDs of metabolites within 4 steps from methionine
KEGG_D5=unique(model.metKeggID(NeighborPos{6}));        %KEGG IDs of metabolites within 5 steps from methionine
Met_D3=unique(K2M_Met(ismember(K2M_KEGG,KEGG_D3)));     %Names of metabolites within 3 steps from methionine
Met_D4=unique(K2M_Met(ismember(K2M_KEGG,KEGG_D4)));     %Names of metabolites within 4 steps from methionine
Met_D5=unique(K2M_Met(ismember(K2M_KEGG,KEGG_D5)));     %Names of metabolites within 5 steps from methionine

%-------------------------------------------------------------------
% Write the metabolite lists to files
%-------------------------------------------------------------------
fileList={'MS_Met3.txt','MS_Met4.txt','MS_Met5.txt'};
MetList={Met_D3,Met_D4,Met_D5};
for i=1:3
    file=fopen(fileList{i},'w');
    for j=1:length(MetList{i})
        fprintf(file,'%s\n',MetList{i}{j});
    end
    fclose(file);
end

%-------------------------------------------------------------------
% Function for identifying neighbors of a metabolite in the genome-
% scale metabolic network Recon 2
%-------------------------------------------------------------------
% INPUT arguments:
% xi - the index for the target metabolite
% S - stoichiometric matrix of the metabolic network
%-------------------------------------------------------------------
function pos=findOtherEnd(xi,S)
    FluxPos=find(S(xi,:)~=0);
    SubMat=S(:,FluxPos);
    SubMat=SubMat.*repmat(S(xi,FluxPos),size(SubMat,1),1);
    pos=find(min(SubMat')<0);
end