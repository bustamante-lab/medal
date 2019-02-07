################################################################################
## Medal algorithm Described in: Medication Usage Similarity in Patients with ##
## Acute-onset Neuropsychiatric Syndrome (Pineda, et. al.)                    ##
##                                                                            ##
## Implemented by Armin Pourshafeie                                           ##
##                                                                            ##
##                                                                            ##
################################################################################


import numpy as np
from numba import  jit,  njit, prange
import numba as nb
import pandas as pd
import time
import tqdm
import sys
import os


START= 1
END  = -1
CONT = 2
DASH = -2

def dataQC(dataName):
    """
    Reads the data into a dataframe and imputes the missing values when possible.
    Removes the rows that cannot be imputed. Finally, it merges rows with overlapping
    times.

    Args:
        dataName: Path to the file with the usage data.

    Returns:
        dataframe with the cleaned up usage data."""

    events = pd.read_csv('medsEvents.csv')
    data = events.loc[:, ["id", "medication", "start", "end"]]
    del events
    # Deal with missing data
    data.loc[data.start < 1, "start"] = 1
    data.loc[data.end < 1, "end"] = 1
    data.loc[np.isnan(data.end), "end"] = data.start + 1
    data.loc[np.isnan(data.start), "start"] = data.end - 1
    data = data.loc[np.logical_not(np.logical_or(np.isnan(data.start), np.isnan(data.end))),:]
    #data.end.loc[data.start == data.end] = data.start[data.start == data.end] + 1

    data.loc[:,"start"] = data.start.astype(int)
    data.loc[:, "end"] = data.end.astype(int)

    # Deal with name changes
    data.loc[data.medication=="Prednisone burst", "medication"] = "Prednisone"
    data.loc[data.medication=="Maintenance prednisone", "medication"] = "Prednisone"

    # combine when time ranges match
    data = range_merger(data, 'start', 'end')
    data.loc[:,"start"] = data.start.astype(int)
    data.loc[:, "end"] = data.end.astype(int)

    return data


#data = pd.DataFrame(data={'id': [1, 1, 1, 1, 1, 2, 2, 2, 2],
#    'medication': ["clindamycin", "clindamycin", "amoxicillin", "amoxicillin", "amoxicillin",
#        "clindamycin", "clindamycin", "amoxicillin", "amoxicillin"],
#    'start': [1, 6, 2, 5, 8, 3, 8, 4, 7],
#    'end'  : [3, 8, 3, 6, 9, 5, 10, 5, 9]
#    }
#    )


#### Bunch of helper function definitions.... can be moved to another file for readability...

def range_merger(df, range_start, range_end):
    """
    Reads the data into a dataframe and imputes the missing values when possible.
    Removes the rows that cannot be imputed. Finally, it merges rows with overlapping
    times.

    Args:
        df         : Usage dataframe.
        range_start: Column name for the start time.
        range_end  : Column name for the end time.

    Returns:
        Dataframe with overlapping times merged.
    """

    names = df.columns
    all_but_range_cols = list(set(names) - set([range_end, range_start]))
    all_but_range_cols = [name for name in names if name in all_but_range_cols]
    df = df.sort_values(by=all_but_range_cols+[range_start]).reset_index(drop=True)
    start = df.loc[0, range_start]
    end   = df.loc[0, range_end]
    merged = []
    prev = df.loc[0,all_but_range_cols]
    for ind, row in df.loc[1:].iterrows():
        if (row.loc[all_but_range_cols] != prev).any():
            merged.append(prev.tolist() + [str(start)] + [str(end)])
            prev = row.loc[all_but_range_cols]
            start = row.loc[range_start]
            end   = row.loc[range_end]
        else: # check range
            if row.loc[range_start] <= end:
                end = max(end, row.loc[range_end])
            else:
                merged.append(prev.tolist() + [str(start)] + [str(end)])
                prev = row.loc[all_but_range_cols]
                start = row.loc[range_start]
                end   = row.loc[range_end]

    # don't forget the last one!
    merged.append(prev.tolist() + [str(start)] + [str(end)])
    df = pd.DataFrame(merged)
    df.columns = all_but_range_cols + [range_start, range_end]
    return df



def _getSequence(pat):
    """
    Encodes a patient usage data into sequence vectors.
    Relies on START, CONT, END values defined globally.

    Args:
        pat: Usage dataframe for a single patient.

    Returns:
        List of unique medications.
        List of sequences (list of numpy arrays).
        List of start times.
        List of end times.
    """
    meds = np.unique(pat.medication)
    sequences = []
    startTimes = []
    endTimes = []
    for medication in meds:
        treated_by_med = pat.loc[pat.medication==medication, :]
        endTime = np.max(treated_by_med.end)
        startTime = np.min(treated_by_med.start)
        startTimes.append(startTime)
        endTimes.append(endTime)
        seq = np.zeros(( endTime - startTime +1), dtype=np.int8)
        base = startTime
        l = len(seq) - 2
        for i_ind, j_ind in zip(treated_by_med.start, treated_by_med.end):
            i = i_ind - base
            j = j_ind - base
            j = min(j, l)
            seq[i] = START
            seq[j+1] = END
            seq[(i+1):j+1] = CONT
        sequences.append(seq)
        seq[-1] = CONT # I disagree with this but...
        if i_ind == j_ind:
            seq[i_ind - base] = START
    return(meds, sequences, startTimes, endTimes)

def getSequences(data):
    """
    Encodes the usage data into sequence vectors.
    Reliies on START, CONT, END values defined globally.

    Args:
        data: Dataframe with all patient's usage data.

    Returns:
        List with tuples of medicatin, sequence, start time and end time for
        each patient
    """

    patients = sorted(np.unique(data.id))
    #npatients = len(patients)
    patientInfo = []
    for patient in patients:
        patientInfo.append(_getSequence(data.loc[data.id ==patient,:]))

    return patientInfo

@njit(nb.int32[:,:](nb.int8[:], nb.int8[:], nb.int32[:,:]), nogil=True,  parallel=True, cache=True)
def paths(seq1, seq2, mat):
    """
    Uses the medal algorithm to populate a cost matrix for different paths.

    Args:
        seq1: Patient1's usage sequence.
        seq2: Patient2's usage sequence.
        mat : int32[:,:] array of size seq1 + 1 x seq1 + 1

    Returns:
        Matrix of cost of various alignments.
    """
    mat[0,:] = 0
    mat[:,0] = 0
    for i in prange(len(seq1)):
        s1 = seq1[i]
        for j in range(len(seq2)):
            s2 = seq2[j]
            #neighbors = [ mat[i,j+1], mat[i,j], mat[i+1,j]] # t, d, l
            if s1 == s2: # Same case
                # WORK AROUND due to numba issues with min(arr)
                val = min(mat[i,j+1], mat[i,j], mat[i+1,j])
            else: # Scenario Changes
                if (s1+s2) % 2 == 0 : # Opposing Scenarios
                    val = max(mat[i,j+1], mat[i,j], mat[i+1,j]) + 1
                elif s1 + s2 == 1 and s1*s2 == 0: # Gap cases
                    val = min(mat[i,j+1], mat[i,j], mat[i+1,j]) + 1
                else:
                    val = max(mat[i,j+1], mat[i,j], mat[i+1,j])
            mat[i+1, j+1] = val
    return (mat)


@jit(nopython=True, nogil=True, cache=True)
def align(seq1, seq2, mat):
    """
    Traceback to determine alignment.

    Args:
        seq1: Expanded usage sequence for patient 1.
        seq2: Expanded usage sequence for patient 2.
        mat : medal alignment path costs.

    Returns:
        alrignment for sequences 1,2 and the cost of the alignment.

    """
    i, j = len(seq1) - 1, len(seq2) - 1
    boolseq1 = seq1 > 0
    boolseq2 = seq2 > 0
    dist   = 0
    align1 = []
    align2 = []
    while (i >= 0 and j >= 0):
        if(boolseq1[i] == boolseq2[j]):
            align1.append(seq1[i])
            align2.append(seq2[j])
            i -= 1
            j -= 1
        elif mat[i, j-1] < mat[i-1, j]:
            align1.append(DASH)
            dist += 1
            align2.append(seq2[j])
            j -= 1
        else:
            align1.append(seq1[i])
            align2.append(DASH)
            dist += 1
            i -= 1
    if i >= 0:
        dist += i+1
        #align1 += seq1[i::-1].tolist()  # append the leftovers
        align1 += [item for item in seq1[i::-1]]  # append the leftovers
        align2 += [DASH for k in range(i+1)]
    elif j >= 0:
        dist += j+1
        #align2 += seq2[j::-1].tolist()  # append the leftovers
        align2 += [item for item in seq2[j::-1]]  # append the leftovers
        align1 += [DASH for k in range(j+1)]
    align1.reverse()
    align2.reverse()
    dist /= 2
    return(align1, align2, dist)



def medalDistance(patientInfo):
    """
    Computes the insertion edit distance between every pair of individuals.

    Args:
        patientInfo: List of meds, sequences, start and end times as produced by getSequence

    Returns:
        Pairwise distance matrix.

    """

    npatients = len(patientInfo)
    total_count = (npatients * (npatients - 1))/(2)
    update_step = int(total_count/100)
    distMat = np.zeros((npatients, npatients))
    pbar = tqdm.tqdm(total = 100)
    counter = 0
    for i in range(npatients):
        meds1, sequences1, startTimes1, endTimes1 = patientInfo[i]
        for j in range(0, i):
            dist, size = 0, 0
            meds2, sequences2, startTimes2, endTimes2 = patientInfo[j]
            inMed2 = meds2.tolist()
            for medInd, med in enumerate(meds1):
                if med not in meds2:
                    days  = endTimes1[medInd] - startTimes1[medInd] + 1
                    dist += days * days
                    size += days
                else:  # It's in both
                    medInd2      = np.where(meds2 == med)[0][0]
                    startTime    = min(startTimes1[medInd], startTimes2[medInd2])
                    endTime      = max(endTimes1[medInd], endTimes2[medInd2])
                    diff         = endTime - startTime + 1
                    expandedSeq1 = np.zeros(diff, dtype=np.int8)
                    expandedSeq2 = np.zeros(diff, dtype=np.int8)
                    expandedSeq1[startTimes1[medInd]-startTime:endTimes1[medInd]-startTime + 1] = sequences1[medInd]
                    expandedSeq2[startTimes2[medInd2]-startTime:endTimes2[medInd2]-startTime + 1] = sequences2[medInd2]
                    if endTimes1[medInd] < endTime:
                        expandedSeq1[endTimes1[medInd] - startTime + 1] = END
                    elif endTimes2[medInd2] < endTime:
                        expandedSeq2[endTimes2[medInd2] - startTime + 1] = END
                    mat = np.empty((diff + 1, diff + 1), dtype=np.int32)
                    mat = paths(expandedSeq1, expandedSeq2, mat)
                    al1, al2, d = align(expandedSeq1, expandedSeq2, mat)
                    size       += len(al1)
                    dist       += d * len(al1)
                    inMed2.remove(med)
            for med in inMed2:
                medInd2 = np.where(meds2 == med)[0][0]
                days    = endTimes2[medInd2] - startTimes2[medInd2] + 1
                dist   += days * days
                size   += days
            distMat[i,j] = dist/float(size)
            #print(dist/float(size))
            counter += 1
            if counter == update_step:
                pbar.update(1)
                counter = 0
    distMat += distMat.T
    return distMat


def medal(usageFile, idFileName="patientID.txt", distMatName="distance_mat.txt"):
    """
    Performs the medal alignment algorithms.

    Args:
        usageFile  : Drug usage file.
        idFileName : File path for the patient ids.
        distMatName: File path for the pairwise distance matrix.

    Returns: None. Writes the patient ids and distance matrix from the medal algorithm.
    """

    t = time.time()
    data = dataQC(usageFile)
    patientInfo = getSequences(data)
    tPreprocess = time.time() - t

    patients = np.unique(data.id)
    del data
    np.savetxt("patientID.txt", patients, fmt='%i')
    #npatients = len(patients)

    tAlign = time.time()
    distMat = medalDistance(patientInfo)
    tAlign = time.time() - tAlign
    np.savetxt("distance_mat.txt", distMat)
    print("Preprocessing time: {}".format( tPreprocess ))
    print("Alignment time: {}".format( tAlign ))
    print("Total time: {}".format( time.time()-t ))

if __name__=="__main__":
    args = sys.argv
    usage = """USAGE: python pymedal datalocation\n """
    if len(args) < 2 or len(args) > 2:
        print (usage)
        print ("""pymedal only takes a single argument""")
    dataFile = args[1]
    if not os.path.isfile(dataFile):
        print (usage)
        print ("{} is not a file".format(dataFile))

    medal(dataFile)
