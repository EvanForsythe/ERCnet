import numpy as np
import pandas as pd

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
