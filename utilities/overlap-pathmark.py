#!/usr/bin/env python
"""overlapPATHMARK.py: 

Usage:
  overlapPATHMARK.py [options] sifFile featureFile

Options:
  -q            run quietly
  -i str        IPL file
"""
## Written By: Sam Ng
## Last Updated: ##/##/####
import getopt, os, re, sys
import mData, mPathway

verbose = True
globalList = "/projects/sysbio/users/TCGA/GBM_MANUSCRIPT/CYTOSCAPE/EXPR_CLUSTER/output.gmt"
#globalList = "/hive/users/sng/pathmark/tcgaKIRC.manuscript/ALTERED/lists_t.tab"
globalPathway = "/projects/sysbio/users/TCGA/GBM_MANUSCRIPT/CYTOSCAPE/EXPR_CLUSTER//pid_110725_pathway.tab"

def usage(code = 0):
    print __doc__
    if code != None:
        sys.exit(code)

def log(msg, die = False):
    if verbose:
        sys.stderr.write(msg)
    if die:
        sys.exit(1)

def syscmd(cmd):
    log("running:\n\t"+cmd+"\n")
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % (exitstatus)
        sys.exit(10)
    log("... done\n")

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "qi:")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 2:
        log("ERROR: incorrect number of arguments", die = True)
    
    sifFile = args[0]
    featureFile = args[1]    
    
    iplMatrix = None
    global verbose, globalList, globalPathway
    for o, a in opts:
        if o == "-q":
            verbose = False
        elif o == "-i":
            iplMatrix = a
    
    ## build sourceList
    typeFile = "/".join(re.split("/", sifFile)[:-2]+["TYPE.NA"])
    sourceList = re.split("/", sifFile)[-2]+".list_t"

    h = mData.rList(featureFile)
    (n, i) = mPathway.rSIF(sifFile, typef = typeFile)
    (gn, gi) = mPathway.rPathway(globalPathway)
    p = mPathway.Pathway(n, i)
    s = mPathway.sortConnected(p)
    f = open(sourceList, "w")
    c = 1
    for i in s:
        u = list(set(i) & set(h))
        if len(u) >= 5:
            f.write("component_%s\t%s\n" % (c, "\t".join(u)))
            c += 1
            break
    f.write("component_all\t%s\n" % ("\t".join(list(set(n.keys()) & set(h)))))
    f.write("all\t%s\n" % ("\t".join(list(set(gn.keys()) & set(h)))))
    f.close()
    
    ## build overlapList
    overlapTotal = 0
    overlapList = "overlap.list_t"
    f = open(globalList, "r")
    o = open(overlapList, "w")
    for line in f:
        line = line.rstrip("\n\r")
        pline = re.split("\t", line)
        pathwayName = pline[0]
        pathwayFeatures = pline[1:]
        o.write("%s\t%s\n" % (pathwayName, "\t".join(list(set(pathwayFeatures) & set(h)))))
	overlapTotal += 1
    overlapTotal = overlapTotal*2
    f.close()
    o.close()
    syscmd("sets_overlap.pl -p 0.05 -I -o \"stats\" %s %s > %s" % (sourceList, overlapList, re.split("/", sifFile)[-2]+".stats.output"))
    syscmd("sets_overlap.pl -p 0.05 -I -o \"members\" %s %s > %s" % (sourceList, overlapList, re.split("/", sifFile)[-2]+".members.output"))
    
    ## read stats output
    statsMap = {}
    f = open(re.split("/", sifFile)[-2]+".stats.output", "r")
    for line in f:
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            pathwayName = line.lstrip(">")
            minVal = 1.0
            pVal = None
        else:
            pline = re.split(",", line)
            if (10**float(pline[1]) < minVal) and (float(pline[2]) >= 5):
                minVal = 10**float(pline[1])
                statsMap[pathwayName] = (pline[0], minVal*float(overlapTotal))
    f.close()
    
    ## read members output
    membersMap = {}
    totalMembers = 0
    proportionMap = {}
    f = open(re.split("/", sifFile)[-2]+".members.output", "r")
    for line in f:
        line = line.rstrip("\n\r")
        pline = re.split("\t", line)
        if pline[1] in statsMap:
            if pline[1] not in membersMap:
                membersMap[pline[1]] = []
                totalMembers = 0
            if pline[2] != "NaN":
                membersMap[pline[1]].append(pline[0])
            totalMembers += 1
            proportionMap[pline[1]] = "%s/%s" % (len(membersMap[pline[1]]), totalMembers)
    f.close()
    
    ## output summary
    f = open(re.split("/", sifFile)[-2]+".summary.tab", "w")
    features = statsMap.keys()
    features.sort(lambda x, y: cmp(statsMap[x][1],statsMap[y][1]))
    for feature in features:
        f.write("%s\t%s\t%s\t%s\t%s\n" % (feature, statsMap[feature][1], ",".join(membersMap[feature]), proportionMap[feature], statsMap[feature][0]))
    f.close()
    
    ## summarize pathways
    if iplMatrix is not None:
        pathwayVals = {}
        pathwayScores = {}
        (iplData, iplSamples, iplFeatures) = mData.rCRSData(iplMatrix, retFeatures = True)
        f = open(re.split("/", sifFile)[-2]+".members.output", "r")
        for line in f:
            line = line.rstrip("\n\r")
            pline = re.split("\t", line)
            if pline[1] in statsMap:
                if pline[1] not in pathwayVals:
                    pathwayVals[pline[1]] = {}
                    for sample in iplSamples:
                        pathwayVals[pline[1]][sample] = []
                if pline[0] in iplFeatures:
                    for sample in iplSamples:
                        pathwayVals[pline[1]][sample].append(abs(float(iplData[sample][pline[0]])))
        f.close()
        for pathwayName in pathwayVals.keys():
            pathwayScores[pathwayName] = {}
            for sample in iplSamples:
                pathwayScores[pathwayName][sample] = mData.mean(pathwayVals[pathwayName][sample])
        mData.wCRSData(re.split("/", sifFile)[-2]+".score", pathwayScores)

if __name__ == "__main__":
    main(sys.argv[1:])
