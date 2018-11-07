import sys,os
cond = False
for line in open(sys.argv[1]):
    if not line.startswith('#'):
        line = line.split('\t')
        if cond:
            print line
            sys.exit(1)
        if line[1] == '21530509':
            print line
            cond = True
