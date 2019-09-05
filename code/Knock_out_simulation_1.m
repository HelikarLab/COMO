% Cobra toolbox
% Create the model
% data =tdfread('input_COBRA_list.txt')
% test_tissue_model = createTissueSpecificModel(model,data)

%import known dt genes (Entrez)
%DT=tdfread('Targets_for_inhibitorsEntrez.txt')
%DT_genes=num2str(DT.Target_ofInhibitos_ENTREZ)
%DT_genes=cellstr(DT_genes)

fid = fopen('Th1_inhibitors_Entrez.txt'); % change filename
DT_genes = textscan(fid,'%s','Delimiter','\n');
DT_genes = DT_genes{1,1};
DT_model=intersect(model.genes,DT_genes);

% reactions of drug target genes
%reaction indices and gene indices from model.rxnGeneMat
[rxnInd, geneInd] = find(model.rxnGeneMat);
% geneInd to gene
geneInd2genes=model.genes(geneInd);
%gtoKD=intersect(DT_model,geneInd2genes)

% simulate WT and gene deletion by perturbing DT_model
WT_sol=optimizeCbModel(model)
WT_sol=WT_sol.x
[grRatio,grRateKO,grRateWT,hasEffect,delRxns] = singleGeneDeletion(model,'MOMA',DT_model)
hasEffect_DT=find(hasEffect==1)
hasEffect_DTgenes=DT_model(hasEffect_DT)
[grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model,'MOMA',hasEffect_DTgenes)

% flux ratio matrix 

fluxSolutionRatios=[];
for i= 1:size(fluxSolution,2)
FSratios=fluxSolution(:,i)./WT_sol;
fluxSolutionRatios(:,i)=FSratios
end

% Read files of up- and down-regulated genes from DAG_genes.txt
%%%%DAG_dis=tdfread('UP_DOWN_DAG.txt')
% Search gene indices for up down genes in diseases
%%DAG_dis_genes=DAG_dis.Gene
%%%DAG_dis_genes= cellstr(DAG_dis_genes)
%%%DAG_dis_genes=strtrim(DAG_dis_genes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('RA_DOWN.txt'); % change filename, run this for both up and down regulated genes
DAG_dis_genes = textscan(fid,'%s','Delimiter','\n');
DAG_dis_genes = DAG_dis_genes{1,1};

% intersect with geneInd2gene to obtain DAG genes within the model
DAG_dis_met_genes=intersect(DAG_dis_genes,geneInd2genes)

% Obtain reactions associated with DAG_dis_met_genes
% first obtain indices from geneIndtogenes
DAG_dis_met_genesInd=find(ismember(geneInd2genes,DAG_dis_met_genes))
DAG_dis_met_rxnInd = rxnInd(DAG_dis_met_genesInd)

% obtain DAG reactions flux ratios
DAG_rxn_fluxRatio = fluxSolutionRatios(DAG_dis_met_rxnInd,:)

%
combined_output={geneInd2genes(DAG_dis_met_genesInd), DAG_dis_met_rxnInd, DAG_rxn_fluxRatio}

% create network of genes based on the change 

gene_mat_out=[];
for i = 1:length(hasEffect_DTgenes)
	FR_i=DAG_rxn_fluxRatio(:,i)
	naN_ind= find(isnan(DAG_rxn_fluxRatio(:,i)));
	FR_i([naN_ind])=[];
	Gene_i=combined_output{1};
	Gene_i([naN_ind])=[];
	Rind_i=combined_output{2};
	Rind_i([naN_ind])=[];

	P_gene=repelem(hasEffect_DTgenes(i),length(Gene_i));
	P_gene=P_gene';
	pegene_mat=[P_gene Gene_i num2cell(FR_i)];
	gene_mat_out=[gene_mat_out; pegene_mat];
end
% Change_file_name here
T = cell2table(gene_mat_out);
writetable(T,'Gene_Pairs_Inhi_Fratio_DOWN.txt');


% get Flux distribution for Rxns of genes associated with Up-regulated genes- score ((total number of flux down/total numebr of fluxes)-total number of upregulated fluxes/total number of fluxes)
% get the flux distribution for reactions associated with down-regulated genes- score 1 if it upregulates the down-regulated gene (total number of down-regulated fluxes/total fluxes-otal number of down-regulated fluxes/total flux)
% get the growth for each knock-out phenotype - score 1 if it not drastically changes the growth 
% consider toxicity - score 1 if it is not changes housekeeping genes
% if it not changes the essential genes - 
% sort scores 
% pick highest one
%


