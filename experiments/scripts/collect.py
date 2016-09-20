import re
import sys
import os

# loop through all files in the directory

# read one file
directory = sys.argv[1]
for fname in os.listdir(directory):
    f = open(directory + "/" + fname)

    durations = [re.findall('^Duration \(s\): .*', line) for line in f]
    durations = [line for line in durations if line]
    output = []
    for e in durations:
        m = re.search("Duration \(s\): ", e[0])
        output.append(e[0][:m.start()] + e[0][m.end():])

    o_file = open('./outputs/' + fname, 'w')
    for o in output:
        o_file.write(o + "\n")
