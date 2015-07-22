import os, sys

# rememeber that for stablematching circuits we must double the inputs size
# because the first step is combining XOR-shared values
# and as of 5/27/15, the latest stablematching circuit requires 32 bits per input before paring down

len = sys.argv[1]

for i in xrange(0,int(len)/4):
    sys.stdout.write('F')
print ""
for i in xrange(0,int(len)/4):
    sys.stdout.write('0')
print ""
