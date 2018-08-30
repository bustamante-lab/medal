import numpy as np
import pdb
from numba import autojit, jit
import pandas as pd
import time

START= 1
END  = -1
CONT = 2

events = pd.read_csv('medsEvents.csv')
meds   = pd.read_csv('medsDictonary.csv')
data = events.loc[:, ["id", "medication", "start", "end"]]
del events
data = pd.DataFrame(data={'id': [1, 1, 1, 1, 1, 2, 2, 2, 2],
    'medication': ["clindamycin", "clindamycin", "amoxicillin", "amoxicillin", "amoxicillin",
        "clindamycin", "clindamycin", "amoxicillin", "amoxicillin"],
    'start': [1, 6, 2, 5, 8, 3, 8, 4, 7],
    'end'  : [3, 8, 3, 6, 9, 5, 10, 5, 9]
    }
    )

### QC
patients = sorted(np.unique( data.id))
npatients = len(patients)

# Deal with missing data
data.loc[data.start < 1, "start"] = 1
data.loc[data.end < 1, "end"] = 1
data.loc[np.isnan(data.end), "end"] = data.start + 1
data.loc[np.isnan(data.start), "start"] = data.end - 1
data = data.loc[np.logical_not(np.logical_or(np.isnan(data.start), np.isnan(data.end))),:]

data.loc[:,"start"] = data.start.astype(int)
data.loc[:, "end"] = data.end.astype(int)

# Deal with name changes
data.loc[data.medication=="Prednisone burst", "medication"] = "Prednisone"
data.loc[data.medication=="Maintenance prednisone", "medication"] = "Prednisone"

#### Bunch of helper function definitions.... can be moved to another file for readability...



def range_merger(df, range_start, range_end):
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

data = range_merger(data, 'start', 'end')
data.loc[:,"start"] = data.start.astype(int)
data.loc[:, "end"] = data.end.astype(int)


def getSequence(pat):
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
        seq = np.zeros(( endTime - startTime +1))
        base = startTime
        l = len(seq) - 2
        for i , j in zip(treated_by_med.start, treated_by_med.end):
            i -= base
            j -= base
            j = min(j, l)
            seq[i] = START
            seq[j+1] = END
            seq[(i+1):j+1] = CONT
        sequences.append(seq)
        seq[-1] = CONT
    return(meds, sequences, startTimes, endTimes)

t = time.time()
medict = {}
patientInfo = []
for patient in patients:
    patientInfo.append(getSequence(data.loc[data.id ==patient,:]))
print("{}".format(time.time() - t))


@jit
def paths(seq1, seq2, mat, dirmat):
    for i, s1 in enumerate(seq1):
        for j, s2 in enumerate(seq2):
            val = 0
            dirs = ""
            neighbors = [ mat[i,j+1], mat[i,j], mat[i+1,j]] # t, d, l
            if s1 == s2: # Same case
                val = min(neighbors)
                dirs = "d"
            else: # Scenario Changes
                if (s1+s2) % 2 == 0 : # Opposing Scenarios
                    mInd = np.argmax(neighbors)
                    val = neighbors[mInd] + 1
                elif s1 + s2 == 1 and s1*s2 == 0: # Gap cases
                    mInd = np.argmin(neighbors)
                    val = neighbors[mInd] + 1
                else:
                    mInd = np.argmax(neighbors)
                    val = neighbors[mInd]
                dirs = ["l", "d", "t"][mInd]

            mat[i+1, j+1] = val
            dirmat[i+1, j+1] = dirs
    return (mat, dirmat)

@jit
def medalDistance(seq1, seq2):
    l1, l2 = len(seq1)+1, len(seq2)+1
    mat = np.zeros((l1, l2))
    dirmat = np.chararray((l1, l2))
    mat, dirmat = paths(seq1, seq2, mat, dirmat)
    return mat, dirmat


distMat = np.zeros((npatients, npatients))

for i in range(0,npatients):
    meds1, sequences1, startTimes1, endTimes1 = patientInfo[i]
    for j in range(0, i):
        dist, size = 0, 0
        meds2, sequences2, startTimes2, endTimes2 = patientInfo[j]
        inMed2 = meds2.tolist()
        for medInd, med in enumerate(meds1):
            if med not in meds2:
                days = endTimes1[medInd] - startTimes1[medInd] + 1
                dist += days
                size += days
            else:  # It's in both
                medInd2 = np.where(meds2 == med)[0][0]
                startTime = min(startTimes1[medInd], startTimes2[medInd2])
                endTime = max(endTimes1[medInd], endTimes2[medInd2])
                diff = endTime - startTime + 1
                expandedSeq1 = np.zeros(diff)
                expandedSeq2 = np.zeros(diff)
                expandedSeq1[startTimes1[medInd]-startTime:endTimes1[medInd]-startTime + 1] = sequences1[medInd]
                expandedSeq2[startTimes2[medInd2]-startTime:endTimes2[medInd2]-startTime + 1] = sequences2[medInd2]
                print(i, j)
                dist, size = medalDistance(expandedSeq1, expandedSeq2)
                inMed2.remove(med)
            pdb.set_trace()
            distMat[i,j] = dist
        for med in inMed2:
            medInd2 = np.where(meds2 == med)[0][0]
            days = endTimes2[medInd2] - startTimes2[medInd2] + 1
            dist += days
            size += days
        distMat[i,j] = dist/float(size)
distMat += distMat.T





