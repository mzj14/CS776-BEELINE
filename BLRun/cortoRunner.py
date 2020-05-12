import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for CORTO.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("CORTO").exists():
        print("Input folder for CORTO does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("CORTO").mkdir(exist_ok = False)
        
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    GeneNames = ExpressionData.index

    geneNum = len(GeneNames)

    for idx in range(geneNum):
        exprPath = "CORTO/ExpressionData" + str(idx) + ".csv"
        SingleCopy = ExpressionData.iloc[[idx]]
        singleName = SingleCopy.index[0]
        SingleCopy.rename(index = {singleName: "UnknownGene"}, inplace = True)
        SupplementData = pd.concat([ExpressionData, SingleCopy])
        SupplementData.to_csv(RunnerObj.inputDir.joinpath(exprPath),
                             sep = ',', header  = True, index = True)
        sourcePath = "CORTO/SourceData" + str(idx) + ".txt"
        SourceNames = GeneNames.to_frame().iloc[list(range(idx)) + list(range(idx + 1, geneNum))]
        SourceNames.to_csv(RunnerObj.inputDir.joinpath(sourcePath),
                             sep = ',', header  = False, index = False)
    

def run(RunnerObj):
    '''
    Function to run CORTO algorithm
    '''
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/CORTO/"
    os.makedirs(outDir, exist_ok = True)

    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    GeneNames = ExpressionData.index

    geneNum = len(GeneNames)
    
    nbootstraps = str(RunnerObj.params['nbootstraps'])

    p = str(RunnerObj.params['p'])

    nthreads = str(RunnerObj.params['nthreads'])
    
    for idx in range(geneNum):
        inputExpressionPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/CORTO/ExpressionData" + str(idx) + ".csv"
    
        inputSourcePath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/CORTO/SourceData" + str(idx) + ".txt"

        outPath = 'data/'+ str(outDir) + 'outFile' + str(idx) + '.txt'
        timePath = 'data/'+ str(outDir) + 'time' + str(idx) + '.txt'
        cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/data/ corto:base /bin/sh -c \"time -v -o', timePath, 'Rscript ./runCorto.R',
                         '-e', inputExpressionPath, '-s', inputSourcePath, '-b', nbootstraps, '-p', p, '-t', nthreads, '-o', outPath, '\"'])
        print(cmdToRun)
        os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from CORTO.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/CORTO/"
    
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    geneNames = ExpressionData.index

    # Read output
    OutSubDF = []
    for idx in range(len(geneNames)):
        # Read output
        outFile = 'outFile'+str(idx)+'.txt'
        if not Path(outDir+outFile).exists():
            # Quit if output file does not exist
            print(outDir+outFile+' does not exist, skipping...')
        else:
            OutSubDF.append(pd.read_csv(outDir+outFile, sep = '\t', header = None))
    
    # megre the dataframe by taking the maximum value from each DF
    # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
    outDF = pd.concat(OutSubDF)
    outDF.columns= ['Gene2','Gene1','EdgeWeight']
    outDF.drop(outDF.index[outDF['Gene1'] == 'UnknownGene'], inplace = True)
    outDF = outDF.reindex(columns = ['Gene1','Gene2','EdgeWeight'])

    # Group by rows code is from here:
    # https://stackoverflow.com/questions/53114609/pandas-how-to-remove-duplicate-rows-but-keep-all-rows-with-max-value
    res = outDF[outDF['EdgeWeight'] == outDF.groupby(['Gene1','Gene2'])['EdgeWeight'].transform('max')]
    # Sort values in the dataframe
    finalDF = res.sort_values('EdgeWeight',ascending=False)

    finalDF.to_csv(outDir+'rankedEdges.csv',sep='\t', index = False)
