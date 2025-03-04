import random
import os
import re
from Bio import Phylo
import time
import glob
import numpy as np
import pandas as pd
from datetime import datetime
import warnings

# Function to search for files by prefix and suffix, and log results to a file
def file_counts_log(jobname, out_dir, directory, prefix, suffix, result_file):
    # Search for files in the directory with the given prefix and suffix
    matching_files = [f for f in os.listdir(out_dir+directory) if f.startswith(prefix) and f.endswith(suffix)]
    
    # Count the number of matching files
    num_matching_files = len(matching_files)
    
    #strip the / from the directory string
    directory = directory.rstrip('/')

    # Write the results to the specified result file in append mode
    with open(result_file, "a") as file_handle:
        file_handle.write(f"{jobname}, {directory}, {num_matching_files}\n")

#Function to count the number of lines in output files and log that info
# Function to count lines in a file and log the result
def file_line_count_log(jobname, file_path, result_file):
    try:
        # Open the file and count the number of lines
        with open(file_path, 'r') as f:
            line_count = sum(1 for _ in f)
        
        #Get only the file name (no path or extension)
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        
        # Write the results to the specified result file in append mode
        with open(result_file, "a") as file_handle:
            file_handle.write(f"{jobname}, {file_name}, {line_count-1}\n")
        
        print(f"Logged {line_count} lines for file '{file_path}'.")
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        with open(result_file, "a") as file_handle:
            file_handle.write(f"{file_name}, File not found\n")


def FilterBranchType(data, filterType):
    branchFilter = data.filter(like=filterType).columns
    
    for i in range(len(branchFilter)):
        data.drop(branchFilter[i], axis='columns', inplace=True)


def FilterCorrelationType(data,filterType):
    corrFilter = data.filter(like=filterType).columns

    for i in range(len(corrFilter)):
        data.drop(corrFilter[i], axis='columns', inplace=True)

def FilterSignificance(data, R, P, filterBy):

    oRows = len(data)

    if filterBy == 'spearman':
        print('Filtering results by Spearman values...')
        data = data[:][data['S_R2'] > R]
        data = data[:][data['S_Pval'] < P]
    elif filterBy == 'pearson':
        print('Filtering results by Pearson values...')
        data = data[:][data['P_R2'] > R]
        data = data[:][data['P_Pval'] < P]
    elif filterBy == 'kendall':
        print('Filtering results by Kendalls Tau values...')
        data = data[:][data['K_Tau'] > R]
        data = data[:][data['K_Pval'] < P]
    else:
        print('Filtering results by Spearman, Pearson, and Kendalls Tau values...')
        data = data[:][data['P_R2'] > R]
        data = data[:][data['P_Pval'] < P]

        data = data[:][data['S_R2'] > R]
        data = data[:][data['S_Pval'] < P]

        data = data[:][data['K_Tau'] > R]
        data = data[:][data['K_Pval'] < P]

    eRows = len(data)

    rowsRemoved = oRows - eRows

    print(str(rowsRemoved) + ' data entries removed due to filtering.\n')

    return data

def FilterFDR(data, R, P, filterBy):

    oRows = len(data)
    
    if filterBy == 'spearman':
        print('Filtering results by FDR Spearman values...')
        data = data[:][data['S_FDR_Corrected_Pval'] < P]
        data = data[:][data['S_R2'] > R]
    elif filterBy == 'pearson':
        print('Filtering results by FDR Pearson values...')
        data = data[:][data['P_FDR_Corrected_Pval'] < P]
        data = data[:][data['P_R2'] > R]
    elif filterBy == 'kendall':
        print('Filtering results by FDR Pearson values...')
        data = data[:][data['K_FDR_Corrected_Pval'] < P]
        data = data[:][data['K_Tau'] > R]
    else:
        print('Filtering results by both FDR Spearman, FDR Pearson, and FDR Kendalls Tau values...')
        data = data[:][data['S_FDR_Corrected_Pval'] < P]
        data = data[:][data['P_FDR_Corrected_Pval'] < P]
        data = data[:][data['K_FDR_Corrected_Pval'] < P]
        data = data[:][data['S_R2'] > R]
        data = data[:][data['P_R2'] > R]
        data = data[:][data['K_Tau'] > R]
        
    eRows = len(data)
    rowsRemoved = oRows - eRows

    print(str(rowsRemoved) + ' data entries removed due to filtering.\n')

    return data

def benchmarkTime(fileName, path, stamp, process, timer):
    timer = time.localtime()
    current_time = time.strftime("%a, %d %b %Y %H:%M:%S", timer)
    with open(path + fileName, "a") as bench:
        bench.write("Parallel " + process + " Call " + stamp  + '\t' + str(current_time) + '\n')
    return 0

import pandas as pd
from datetime import datetime
import warnings

def benchmarkProcess(bench_filePath, num_items):
    # Load the benchmark file
    try:
        benchmark_file = pd.read_csv(bench_filePath, delimiter='\t', header=None, names=["Process", "Timestamp"])
    except Exception as e:
        raise ValueError(f"Error reading the benchmark file: {e}")
    
    # Helper function to parse date-time safely
    def parse_datetime(value):
        try:
            return datetime.strptime(value, "%a, %d %b %Y %H:%M:%S")
        except ValueError:
            warnings.warn(f"Invalid date-time format: {value}. Replacing with NA.")
            return None  # Use None to represent missing/invalid data
    
    # Find the most recent "Call start" and "Call end" pair
    start_index, end_index = None, None
    for i in range(len(benchmark_file) - 1, -1, -1):  # Iterate backward
        if "Call end" in benchmark_file.iloc[i]["Process"] and end_index is None:
            end_index = i
        elif "Call start" in benchmark_file.iloc[i]["Process"] and start_index is None:
            start_index = i
            break
    
    if start_index is None or end_index is None:
        warnings.warn("No recent 'Call start' and 'Call end' pair found.")
        return 0  # Exit if no valid pair is found
    
    # Parse the start and end times
    process_name = benchmark_file.iloc[start_index]["Process"].replace("Call start", "").strip()
    t_start = parse_datetime(benchmark_file.iloc[start_index]["Timestamp"])
    t_end = parse_datetime(benchmark_file.iloc[end_index]["Timestamp"])
    
    # Calculate total time and lines per minute
    if t_start is None or t_end is None:
        t_int_minutes, lines_per_min = "NA", "NA"
    else:
        t_total = t_end - t_start
        t_int_minutes = max(t_total.total_seconds() / 60, 0)
        lines_per_min = num_items / t_int_minutes if t_int_minutes > 0 else num_items
        lines_per_min = min(lines_per_min, num_items)
    
    # Append results to the benchmark file
    try:
        with open(bench_filePath, "a") as bench:
            bench.write(f"Process: {process_name}\n")
            bench.write(f"Total time (m)\t{t_int_minutes}\n")
            bench.write(f"Items completed per minute\t{lines_per_min}\n\n")
    except Exception as e:
        raise IOError(f"Error writing to the benchmark file: {e}")
    
    return 0

def CheckFileExists(filePath):
    
    if len(glob.glob(filePath)) > 0:
        return True
    else:
        return False

def CheckFileNonEmpty(filePath):
    return os.path.getsize(filePath) > 0

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


# Function for resolving non-bifurcating trees
def resolve_polytomies(tree, min_branch_length=0.001, max_branch_length=0.01):
    """
    Recursively resolves polytomies (non-bifurcating nodes) by randomly adding internal nodes.
    
    Parameters:
    tree (Phylo.BaseTree.Tree): The input phylogenetic tree.
    min_branch_length (float): Minimum branch length for newly added branches.
    max_branch_length (float): Maximum branch length for newly added branches.
    
    Returns:
    Phylo.BaseTree.Tree: A bifurcating version of the input tree.
    """
    
    def is_polytomy(clade):
        """Check if a clade is a polytomy (has more than two children)."""
        return len(clade.clades) > 2
    
    def resolve_clade(clade):
        """Recursively resolve polytomies in a clade."""
        # If the clade is a polytomy, resolve it
        if is_polytomy(clade):
            while len(clade.clades) > 2:
                # Randomly pick two children to group into a new internal node
                child1, child2 = random.sample(clade.clades, 2)
                clade.clades.remove(child1)
                clade.clades.remove(child2)
                
                # Create a new internal node with the two children
                new_internal = Phylo.BaseTree.Clade(branch_length=random.uniform(min_branch_length, max_branch_length))
                new_internal.clades.append(child1)
                new_internal.clades.append(child2)
                
                # Add the new internal node back to the clade
                clade.clades.append(new_internal)
        
        # Recursively resolve polytomies in each child clade
        for child in clade.clades:
            resolve_clade(child)
    
    # Resolve polytomies starting from the root of the tree
    resolve_clade(tree.root)
    
    return tree

#Due to importing errors, this function has been pulled from scipy.stats base code.
#All credit for development and use of this function goes to the scipy creators. 
#(Hopefully we can delete this and use the scipy library version soon)
def false_discovery_control(ps, *, axis=0, method='bh'):
    """Adjust p-values to control the false discovery rate.

    The false discovery rate (FDR) is the expected proportion of rejected null
    hypotheses that are actually true.
    If the null hypothesis is rejected when the *adjusted* p-value falls below
    a specified level, the false discovery rate is controlled at that level.

    Parameters
    ----------
    ps : 1D array_like
        The p-values to adjust. Elements must be real numbers between 0 and 1.
    axis : int
        The axis along which to perform the adjustment. The adjustment is
        performed independently along each axis-slice. If `axis` is None, `ps`
        is raveled before performing the adjustment.
    method : {'bh', 'by'}
        The false discovery rate control procedure to apply: ``'bh'`` is for
        Benjamini-Hochberg [1]_ (Eq. 1), ``'by'`` is for Benjaminini-Yekutieli
        [2]_ (Theorem 1.3). The latter is more conservative, but it is
        guaranteed to control the FDR even when the p-values are not from
        independent tests.

    Returns
    -------
    ps_adusted : array_like
        The adjusted p-values. If the null hypothesis is rejected where these
        fall below a specified level, the false discovery rate is controlled
        at that level.

    See Also
    --------
    combine_pvalues
    statsmodels.stats.multitest.multipletests

    Notes
    -----
    In multiple hypothesis testing, false discovery control procedures tend to
    offer higher power than familywise error rate control procedures (e.g.
    Bonferroni correction [1]_).

    If the p-values correspond with independent tests (or tests with
    "positive regression dependencies" [2]_), rejecting null hypotheses
    corresponding with Benjamini-Hochberg-adjusted p-values below :math:`q`
    controls the false discovery rate at a level less than or equal to
    :math:`q m_0 / m`, where :math:`m_0` is the number of true null hypotheses
    and :math:`m` is the total number of null hypotheses tested. The same is
    true even for dependent tests when the p-values are adjusted accorded to
    the more conservative Benjaminini-Yekutieli procedure.

    The adjusted p-values produced by this function are comparable to those
    produced by the R function ``p.adjust`` and the statsmodels function
    `statsmodels.stats.multitest.multipletests`. Please consider the latter
    for more advanced methods of multiple comparison correction.

    References
    ----------
    .. [1] Benjamini, Yoav, and Yosef Hochberg. "Controlling the false
           discovery rate: a practical and powerful approach to multiple
           testing." Journal of the Royal statistical society: series B
           (Methodological) 57.1 (1995): 289-300.

    .. [2] Benjamini, Yoav, and Daniel Yekutieli. "The control of the false
           discovery rate in multiple testing under dependency." Annals of
           statistics (2001): 1165-1188.

    .. [3] TileStats. FDR - Benjamini-Hochberg explained - Youtube.
           https://www.youtube.com/watch?v=rZKa4tW2NKs.

    .. [4] Neuhaus, Karl-Ludwig, et al. "Improved thrombolysis in acute
           myocardial infarction with front-loaded administration of alteplase:
           results of the rt-PA-APSAC patency study (TAPS)." Journal of the
           American College of Cardiology 19.5 (1992): 885-891.

    Examples
    --------
    We follow the example from [1]_.

        Thrombolysis with recombinant tissue-type plasminogen activator (rt-PA)
        and anisoylated plasminogen streptokinase activator (APSAC) in
        myocardial infarction has been proved to reduce mortality. [4]_
        investigated the effects of a new front-loaded administration of rt-PA
        versus those obtained with a standard regimen of APSAC, in a randomized
        multicentre trial in 421 patients with acute myocardial infarction.

    There were four families of hypotheses tested in the study, the last of
    which was "cardiac and other events after the start of thrombolitic
    treatment". FDR control may be desired in this family of hypotheses
    because it would not be appropriate to conclude that the front-loaded
    treatment is better if it is merely equivalent to the previous treatment.

    The p-values corresponding with the 15 hypotheses in this family were

    >>> ps = [0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344,
    ...       0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.000]

    If the chosen significance level is 0.05, we may be tempted to reject the
    null hypotheses for the tests corresponding with the first nine p-values,
    as the first nine p-values fall below the chosen significance level.
    However, this would ignore the problem of "multiplicity": if we fail to
    correct for the fact that multiple comparisons are being performed, we
    are more likely to incorrectly reject true null hypotheses.

    One approach to the multiplicity problem is to control the family-wise
    error rate (FWER), that is, the rate at which the null hypothesis is
    rejected when it is actually true. A common procedure of this kind is the
    Bonferroni correction [1]_.  We begin by multiplying the p-values by the
    number of hypotheses tested.

    >>> import numpy as np
    >>> np.array(ps) * len(ps)
    array([1.5000e-03, 6.0000e-03, 2.8500e-02, 1.4250e-01, 3.0150e-01,
           4.1700e-01, 4.4700e-01, 5.1600e-01, 6.8850e-01, 4.8600e+00,
           6.3930e+00, 8.5785e+00, 9.7920e+00, 1.1385e+01, 1.5000e+01])

    To control the FWER at 5%, we reject only the hypotheses corresponding
    with adjusted p-values less than 0.05. In this case, only the hypotheses
    corresponding with the first three p-values can be rejected. According to
    [1]_, these three hypotheses concerned "allergic reaction" and "two
    different aspects of bleeding."

    An alternative approach is to control the false discovery rate: the
    expected fraction of rejected null hypotheses that are actually true. The
    advantage of this approach is that it typically affords greater power: an
    increased rate of rejecting the null hypothesis when it is indeed false. To
    control the false discovery rate at 5%, we apply the Benjamini-Hochberg
    p-value adjustment.

    >>> from scipy import stats
    >>> stats.false_discovery_control(ps)
    array([0.0015    , 0.003     , 0.0095    , 0.035625  , 0.0603    ,
           0.06385714, 0.06385714, 0.0645    , 0.0765    , 0.486     ,
           0.58118182, 0.714875  , 0.75323077, 0.81321429, 1.        ])

    Now, the first *four* adjusted p-values fall below 0.05, so we would reject
    the null hypotheses corresponding with these *four* p-values. Rejection
    of the fourth null hypothesis was particularly important to the original
    study as it led to the conclusion that the new treatment had a
    "substantially lower in-hospital mortality rate."

    """
    # Input Validation and Special Cases
    ps = np.asarray(ps)

    ps_in_range = (np.issubdtype(ps.dtype, np.number)
                   and np.all(ps == np.clip(ps, 0, 1)))
    if not ps_in_range:
        raise ValueError("`ps` must include only numbers between 0 and 1.")

    methods = {'bh', 'by'}
    if method.lower() not in methods:
        raise ValueError(f"Unrecognized `method` '{method}'."
                         f"Method must be one of {methods}.")
    method = method.lower()

    if axis is None:
        axis = 0
        ps = ps.ravel()

    axis = np.asarray(axis)[()]
    if not np.issubdtype(axis.dtype, np.integer) or axis.size != 1:
        raise ValueError("`axis` must be an integer or `None`")

    if ps.size <= 1 or ps.shape[axis] <= 1:
        return ps[()]

    ps = np.moveaxis(ps, axis, -1)
    m = ps.shape[-1]

    # Main Algorithm
    # Equivalent to the ideas of [1] and [2], except that this adjusts the
    # p-values as described in [3]. The results are similar to those produced
    # by R's p.adjust.

    # "Let [ps] be the ordered observed p-values..."
    order = np.argsort(ps, axis=-1)
    ps = np.take_along_axis(ps, order, axis=-1)  # this copies ps

    # Equation 1 of [1] rearranged to reject when p is less than specified q
    i = np.arange(1, m+1)
    ps *= m / i

    # Theorem 1.3 of [2]
    if method == 'by':
        ps *= np.sum(1 / i)

    # accounts for rejecting all null hypotheses i for i < k, where k is
    # defined in Eq. 1 of either [1] or [2]. See [3]. Starting with the index j
    # of the second to last element, we replace element j with element j+1 if
    # the latter is smaller.
    np.minimum.accumulate(ps[..., ::-1], out=ps[..., ::-1], axis=-1)

    # Restore original order of axes and data
    np.put_along_axis(ps, order, values=ps.copy(), axis=-1)
    ps = np.moveaxis(ps, -1, axis)

    return np.clip(ps, 0, 1)
