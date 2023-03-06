import os
import time
import glob
import numpy as np
import pandas as pd
from datetime import datetime


def FilterBranchType(data, filterType):
    branchFilter = data.filter(like=filterType).columns
    
    for i in range(len(branchFilter)):
        data.drop(branchFilter[i], axis='columns', inplace=True)


def FilterCorrelationType(data,filterType):
    corrFilter = data.filter(like=filterType).columns

    for i in range(len(corrFilter)):
        data.drop(corrFilter[i], axis='columns', inplace=True)

def FilterSignificance(data, R, P):

    data = data[:][data['R2'] > R]
    data = data[:][data['Pval'] < P]
    return data

def benchmarkTime(fileName, path, stamp, process, timer):
    timer = time.localtime()
    current_time = time.strftime("%a, %d %b %Y %H:%M:%S", timer)
    with open(path + fileName, "a") as bench:
        bench.write("Parallel " + process + " Call " + stamp  + '\t' + str(current_time) + '\n')
    return 0

def benchmarkProcess(bench_filePath, num_items):
    benchmark_file = pd.read_csv(bench_filePath, delimiter = '\t')

    t_finish = datetime.strptime(str(benchmark_file.iloc[-1,-1]), "%a, %d %b %Y %H:%M:%S")
    t_start = datetime.strptime(str(benchmark_file.iloc[-2,-1]), "%a, %d %b %Y %H:%M:%S")
    t_total = t_finish - t_start
    t_int_minutes = float(t_total.total_seconds() / 60)

    if t_int_minutes > 0:
        lines_per_min = num_items / t_int_minutes
        if lines_per_min > num_items:
            lines_per_min = num_items
    else:
        lines_per_min = num_items

    with open(bench_filePath, "a") as bench:
        bench.write("Total time (m)" + '\t' + str(t_int_minutes) + '\n')
        bench.write("Items completed per minute" + '\t' + str(lines_per_min) + '\n''\n')
    return 0

def CheckFileExists(filePath):
    
    if len(glob.glob(filePath)) > 0:
        return True
    else:
        return False

def GetLatestFileName(filePath, fileName):

    fileNames = os.listdir(filePath)
    filteredFileNames = filter(lambda fn: fileName in fn, fileNames)
    
    filteredFileNames = list(filteredFileNames)
    filteredFileNames.sort()

    return filteredFileNames[-1]

def AppendFileName(fileName):

    index = -1
    start = len(fileName) - 5
    end = len(fileName) - 4
    fileNum = 0

    if (fileName[start:end].isnumeric()):
        for i in range(len(fileName)): 
            index += 1
            start -= 1
            if not(fileName[start:end].isnumeric()):
                break
        
        fileNum = int(fileName[(start + 1):end])
        fileNum += 1
        fileName = fileName[0:(end - (index + 2))] + '_' + str(fileNum) + '.tsv'

    else:
        fileNum = 1
        fileName = fileName[0:end] + '_' + str(fileNum) + '.tsv'      

    return fileName

def CheckAndMakeDir(out_dir, dir_name):

    if not os.path.isdir(out_dir + dir_name):
        os.makedirs(out_dir + dir_name)
        print("Created Folder " + dir_name + '\n')
    
    print("Results will be stored in " + dir_name + '\n\n')

    return
