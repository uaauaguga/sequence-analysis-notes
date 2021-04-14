#!/usr/bin/env python
import re
with open("params/RIBOSUM90-unpaired.txt") as f:
    line = next(f)
    pairs = re.split(r"\s+",line.strip())
    for i,line in enumerate(f):
        fields = re.split(r"\s+",line.strip())
        x = fields[0]
        for y,score in zip(pairs,fields[1:]):
            print(x,y,score,sep="\t")
        
