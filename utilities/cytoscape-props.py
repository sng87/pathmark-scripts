import os, re
import mData

assert(os.path.exists("first.props"))
assert(os.path.exists("include.features"))

features = mData.rList("include.features")

f = open("first.props", "r")
o = open("all.props", "w")
for line in f:
    if line.startswith("nodeCustomGraphics1.default-Node\ Custom\ Graphics\ 1-Discrete\ Mapper.mapping.map"):
        pline = re.split(",", line.rstrip("\r\n"))
        imgCounter = int(pline[1])
        for feature in features:
            reline = re.sub(features[0], feature, pline[0]) + "," + str(imgCounter) + "," + re.sub(features[0], feature, pline[1]) + "," + pline[2]
            imgCounter += 1
    else:
        o.write(line)
f.close()
o.close()