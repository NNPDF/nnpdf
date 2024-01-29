import pandas as pd
import sys
import pathlib

def split_string(test_str):
    substrings = test_str.split()
    substrings = list(filter(lambda s: s.strip(), substrings))
    return substrings

FILE_TO_IMPLEMENT = pathlib.Path(sys.argv[1])
FILE_TO_DUMP = pathlib.Path(sys.argv[2])
TYPE_PROC = sys.argv[3]
NAME = FILE_TO_IMPLEMENT.stem

###############
#READ THE FILE
with open(FILE_TO_IMPLEMENT, "r") as f:
    lines = f.readlines()
Q2s = []
x = []
y = []
cross_sections = []
for line in lines[10:]:
    cleanedline = split_string(line.strip('\n'))
    x.append(float(cleanedline[0]))
    Q2s.append(float(cleanedline[1]))
    y.append(float(cleanedline[2]))
    cross_sections.append(float(cleanedline[3]))
LENGHT = x.__len__()
################
#ADD REAL DATA

TYPE_list = [TYPE_PROC for i in range(LENGHT)]
INDEX_list = [i+1 for i in range(LENGHT)]
STAT_LIST = [0.0 for i in range(LENGHT)]
SYS_LIST_ADD = [0.0 for i in range(LENGHT)]
SYS_LIST_MULT = [0.0 for i in range(LENGHT)]

FINAL_DATAFRAME = pd.DataFrame.from_dict({0:INDEX_list,1:TYPE_list,2:x,3:Q2s,4:y,5:cross_sections,6:STAT_LIST, 7:SYS_LIST_ADD, 8:SYS_LIST_MULT})
with open(FILE_TO_DUMP, 'a') as file:
    file.write(f'{NAME}\t1\t{LENGHT}\n')
    FINAL_DATAFRAME.to_csv(file, sep="\t", header=False, index=False)


