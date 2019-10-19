#!/usr/bin/python3
import sys
import getopt
# sys.path.extend(['E:\\reconstruction', 'E:\\cobamp', 'E:\\reconstruction\\src', 'E:\\cobamp\\src', 'E:/reconstruction'])
import cobra
import framed
import cobamp
import pandas as pd
import numpy as np
import pickle
import scipy
import os
import re
import collections

# for testing the algorithms
from cobamp.wrappers import MatFormatReader
from cobamp.wrappers import COBRAModelObjectReader
from troppo.methods.imat import IMAT
from troppo.methods.gimme import GIMME
from troppo.methods.tINIT import tINIT
from troppo.methods.fastcore import FASTcore
from troppo.reconstruction_properties import FastcoreProperties, tINITProperties, GIMMEProperties, IMATProperties
from troppo.utilities.statistics import normalize, z_score

from project import configs


def float_logical_exp(expressionIn, level=0):
    try:
        loc_r = expressionIn.index(')')
    except:
        if 'and' in expressionIn:
            expressionIn = expressionIn.replace('and', ',')
            expressionIn = 'min{' + expressionIn + '}'
        elif 'or' in expressionIn:
            expressionIn = expressionIn.replace('or', ',')
            expressionIn = 'min{' + expressionIn + '}'
        else:
            expressionIn = expressionIn.replace('[', '')
            expressionIn = expressionIn.replace(']', '')
        return expressionIn
    loc_l = expressionIn[0:loc_r].rindex('(')
    innerstring = expressionIn[loc_l:loc_r+1]
    innerstring = innerstring.replace('(', '[')
    innerstring = innerstring.replace(')', ']')
    if 'and' in innerstring:
        innerstring = innerstring.replace('and', ',')
        innerstring = 'min{'+innerstring+'}'
    elif 'or' in innerstring:
        innerstring = innerstring.replace('or', ',')
        innerstring = 'min{'+innerstring+'}'
    else:
        innerstring = innerstring.replace('[', '')
        innerstring = innerstring.replace(']', '')

    expressionOut = '{}{}{}'.format(expressionIn[0:loc_l], innerstring, expressionIn[loc_r+1:])
    expressionOut = float_logical_exp(expressionOut, level+1)

    return expressionOut


def gene_rule_float(expressionIn):
    gene_reaction_by_rule = float_logical_exp(expressionIn)
    gene_reaction_by_rule = gene_reaction_by_rule.replace('{', '(')
    gene_reaction_by_rule = gene_reaction_by_rule.replace('}', ')')
    return gene_reaction_by_rule


def createTissueSpecificModel(GeneralModelFile, GeneExpressionFile):
    model_cobra = cobra.io.load_matlab_model(GeneralModelFile)

    mat = scipy.io.loadmat(GeneralModelFile)['Teff']
    model = MatFormatReader(mat)
    S = model.S
    lb, ub = model.get_model_bounds(False, True)
    rx_names = model.get_reaction_and_metabolite_ids()[0]

    expressionRxns, expVector = mapExpressionToRxn(model_cobra, GeneExpressionFile)

    idx_objective = rx_names.index('biomass_reaction_Mphage')
    properties = GIMMEProperties(exp_vector=expVector,#np.array(gimme_data['0']),
                                obj_frac=0.9,
                                objectives= [{idx_objective:1}],
                                preprocess=True,
                                flux_threshold=0.9
                                )
    algorithm = GIMME(S, lb.astype(float), ub.astype(float), properties)
    model_GIMME = algorithm.run()

    model_GIMME_final = model_cobra.copy() # this is done since it alters the original model_cobra; this way is to guarantee that a new model is changed instead of the original model
    r_ids = [r.id for r in model_GIMME_final.reactions]
    to_remove_ids = [r_ids[r] for r in np.where(model_GIMME==0)[0]]
    model_GIMME_final.remove_reactions(to_remove_ids,True) # this is to get the ids of the reactions to be removed in the model; True is to remove the pending genes/metabolites that with the removal of the reaction can no longer be connected in the network

    print('1\'s: ' + str(len(np.where(model_GIMME==1)[0])))
    print('2\'s: ' + str(len(np.where(model_GIMME==2)[0])))

    return model_GIMME_final


def splitGeneExpressionFile(GeneExpressionFile):
    expressionData = pd.read_csv(GeneExpressionFile)
    expressionData.rename(columns={'ENTREZ_GENE_ID': 'Gene', 'Express': 'Data'},inplace=True)
    expressionData = expressionData.loc[:, ['Gene', 'Data']]
    singleGeneNames = expressionData[~expressionData.Gene.str.contains('//')].reset_index(drop=True)
    multipleGeneNames = expressionData[expressionData.Gene.str.contains('//')].reset_index(drop=True)
    breaksGeneNames = pd.DataFrame(columns=['Gene', 'Data'])
    print(singleGeneNames.shape)
    print(multipleGeneNames.shape)
    for index,row in multipleGeneNames.iterrows():
        for genename in row['Gene'].split('///'):
            breaksGeneNames = breaksGeneNames.append({'Gene': genename, 'Data': row['Data']},ignore_index=True)
    print(breaksGeneNames.shape)
    GeneExpressions = singleGeneNames.append(breaksGeneNames,ignore_index=True)
    print(GeneExpressions.shape)
    GeneExpressions.set_index('Gene',inplace=True)
    return GeneExpressions


def mapExpressionToRxn(model_cobra, GeneExpressionFile):
    GeneExpressions = splitGeneExpressionFile(GeneExpressionFile)
    # GeneExpressions

    expressionRxns = collections.OrderedDict()
    # expVector = []
    cnt = 0
    for rxn in model_cobra.reactions:
        gene_reaction_rule = rxn.gene_reaction_rule
        gene_ids = re.findall(r'\d+', gene_reaction_rule)
        expressionRxns[rxn.id] = -1.0
        if gene_reaction_rule == '':
            continue
        for id in gene_ids:
            boolval = '-1'
            if id in GeneExpressions.index:
                boolval = '{}'.format(GeneExpressions.at[id, 'Data'])
            gene_reaction_rule = gene_reaction_rule.replace('{}'.format(id), boolval, 1)

        try:
            gene_reaction_by_rule = gene_rule_float(gene_reaction_rule)
            expressionRxns[rxn.id] = eval(gene_reaction_by_rule)
            # expVector.append(expressionRxns[rxn.id])
        except:
            print(gene_reaction_by_rule)
            cnt += 1
    print('Map gene expression to reactions, {} errors.'.format(cnt))
    expVector = np.array(list(expressionRxns.values()),dtype=np.float)
    return expressionRxns, expVector


# testexp = '(( 2977 ) and ( 2983 )) or 4882 or (( 2982 ) and ( 2974 )) or 2986 or 2984 or 4881 or (( 2977 ) and ( 2974 )) or 3000 or (( 2982 ) and ( 2983 ))'

# if not (testexp[0] == '(' and testexp[-1] == ')'):
#     testexp = '('+testexp+')'
# outexp=float_logical_exp(testexp)
#
# outexp = outexp.replace('{','(')
# outexp = outexp.replace('}',')')
# print(eval(outexp))

def main(argv):
    modelfile = 'GeneralModel.mat'
    genefile = 'merged_Th1.csv'
    outputfile = 'Th1_SpecificModel.mat'
    try:
        opts, args = getopt.getopt(argv, "hm:g:o:", ["mfile=", "gfile=", "ofile="])
    except getopt.GetoptError:
        print('python3 create_tissue_specific_model.py -m <modelfile> -g <genefile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 create_tissue_specific_model.py -m <modelfile> -g <genefile> -o <outputfile>')
            sys.exit()
        elif opt in ("-m", "--mfile"):
            modelfile = arg
        elif opt in ("-g", "--gfile"):
            genefile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print('General Model file is "{}"'.format(modelfile))
    print('Gene Expression file is "{}"'.format(genefile))
    print('Output file is "{}"'.format(outputfile))
    #print(configs.rootdir)
    GeneralModelFile = os.path.join(configs.datadir, modelfile)
    GeneExpressionFile = os.path.join(configs.datadir, genefile)
    TissueModel = createTissueSpecificModel(GeneralModelFile, GeneExpressionFile)
    print(TissueModel)
    cobra.io.mat.save_matlab_model(TissueModel, os.path.join(configs.datadir, outputfile))
    # cobra.io.sbml.write_cobra_model_to_sbml_file(TissueModel, os.path.join(configs.datadir,'Th1_SpecificModel.xml'))
    print('Genes: ' + str(len(TissueModel.genes)))
    print('Metabolites: ' + str(len(TissueModel.metabolites)))
    print('Reactions: ' + str(len(TissueModel.reactions)))
    print(TissueModel.objective._get_expression())
    print(TissueModel.optimize())

if __name__ == "__main__":
   main(sys.argv[1:])
