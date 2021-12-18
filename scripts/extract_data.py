from sys import argv
from re import sub

inp = argv[1]
print(argv[1])
with open(inp) as f:
    data = f.readline()
    if not data[0] == '[':
        print("bad file, no list")
        exit(1)
    
    data = data.split('] [')
    data = [sub('\[|\]|,', '', l) for l in data]

outname = sub('\.|\/|out', '', inp)+"_data.txt"
with open(outname, 'w+') as f:
    for l in data:
        f.write(l+'\n')
