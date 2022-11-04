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

def FilterSignificance(data, R, P, filterType, corrType):

    if (filterType == 'BXB'):
        branchMethod = 'R2T'
    else:
        branchMethod = 'BXB'

    if (corrType == 'Pearson'):
        corrMethod = 'Spearman'
    else:
        corrMethod = 'Pearson'

    rFilter = corrMethod + '_R2_' + branchMethod
    pFilter = corrMethod + '_P_' + branchMethod

    data = data[:][data[rFilter] > R]
    data = data[:][data[pFilter] < P]

    return data
