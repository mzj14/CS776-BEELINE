import os
import pandas as pd
from pathlib import Path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for TENET.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("TENET").exists():
        print("Input folder for TENET does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("TENET").mkdir(exist_ok = False)
    
    '''
    if not RunnerObj.inputDir.joinpath("TENET/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        ExpressionData.T.to_csv(RunnerObj.inputDir.joinpath("TENET/ExpressionData.csv"),
                             sep = ',', header  = True, index = True)
    
    if not RunnerObj.inputDir.joinpath("TENET/CellData.txt").exists():
        PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                                     header = 0, index_col = 0)
        PTData.to_csv(RunnerObj.inputDir.joinpath("TENET/CellData.txt"),
                             sep = ',', header  = False, index = False)
        SelectInfo = [1] * PTData.size
        SelectData = pd.DataFrame(SelectInfo)
        SelectData.to_csv(RunnerObj.inputDir.joinpath("TENET/SelectData.txt"),
                             sep = ',', header  = False, index = False)
    '''

    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    
    colNames = PTData.columns
    
    for idx in range(len(colNames)):
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        
        exprName = "TENET/ExpressionData" + str(idx) + ".csv"
        ExpressionData.loc[:,index].T.to_csv(RunnerObj.inputDir.joinpath(exprName),
                             sep = ',', header  = True, index = True)
        
        cellName = "TENET/CellData" + str(idx) + ".txt"
        ptDF = PTData.loc[index, [colName]]
        ptDF.to_csv(RunnerObj.inputDir.joinpath(cellName), 
                sep = ',', header  = False, index = False)
        
        selectName = "TENET/SelectData" + str(idx) + ".txt"
        SelectInfo = [1] * ptDF.size
        SelectData = pd.DataFrame(SelectInfo)
        SelectData.to_csv(RunnerObj.inputDir.joinpath(selectName),
                             sep = ',', header  = False, index = False)

def run(RunnerObj):
    '''
    Function to run TENET algorithm
    '''
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    
    colNames = PTData.columns
    
    number_of_threads = str(RunnerObj.params['number_of_threads'])

    history_length = str(RunnerObj.params['history_length'])

    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/TENET/"
    os.makedirs(outDir, exist_ok = True)


    for idx in range(len(colNames)): 
        inputExpressionPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/TENET/ExpressionData" + str(idx) + ".csv"
    
        inputCellPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/TENET/CellData" + str(idx) + ".txt"

        inputSelectPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/TENET/SelectData" + str(idx) + ".txt"

        outPath = 'data/'+ str(outDir) + 'outFile' + str(idx) + '.txt'
        
        timePath = 'data/'+ str(outDir) + 'time' + str(idx) + '.txt'

        cmdToRun = ' '.join(['docker run --rm -v', str(Path.cwd())+':/runTENET/data/ tenet:base /bin/sh -c \"time -v -o', timePath, './AllinOneRun.sh',
                         inputExpressionPath, number_of_threads, inputCellPath, inputSelectPath, history_length, outPath, '\"'])
        print(cmdToRun)
        os.system(cmdToRun)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from TENET.
    '''
    # Quit if output directory does not exist
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/TENET/"
    
    # Read output
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)
    colNames = PTData.columns
    for idx in range(len(colNames)):
        OutDF = pd.read_csv(outDir+'outFile' + str(idx) + '.txt', sep = '\t', header = 0, index_col = 0)
        outFile = open(outDir + 'outFile' + str(idx) + '.csv','w') 
        # Form ranked edges
        rows = list()
        for jdx in OutDF.index:
            for kdx in OutDF.index:
                outFile.write('\t'.join([kdx,jdx,str(OutDF[jdx][kdx])])+'\n')
        outFile.close()

    OutSubDF = [0]*len(colNames)
    for idx in range(len(colNames)):
        # Read output
        outFile = 'outFile'+str(idx)+'.csv'
        if not Path(outDir+outFile).exists():
            # Quit if output file does not exist

            print(outDir+outFile+' does not exist, skipping...')
            return
        OutSubDF[idx] = pd.read_csv(outDir+outFile, sep = '\t', header = None)

    # megre the dataframe by taking the maximum value from each DF
    # From here: https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
    outDF = pd.concat(OutSubDF)
    outDF.columns= ['Gene1','Gene2','EdgeWeight']
    # Group by rows code is from here:
    # https://stackoverflow.com/questions/53114609/pandas-how-to-remove-duplicate-rows-but-keep-all-rows-with-max-value
    res = outDF[outDF['EdgeWeight'] == outDF.groupby(['Gene1','Gene2'])['EdgeWeight'].transform('max')]
    # Sort values in the dataframe
    finalDF = res.sort_values('EdgeWeight',ascending=False)

    finalDF.to_csv(outDir+'rankedEdges.csv',sep='\t', index = False)
    
