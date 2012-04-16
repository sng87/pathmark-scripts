#!/usr/bin/env python
"""statistics-subnets.py:

Usage:
  statistics-subnets.py [options] < nets

Options:
  -c str        counts stat file
  -q            run quietly
"""
import os, os.path, sys, getopt, re, math
import mData, mPathway, mCalculate

verbose = True
outStats = False

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg):
    if (verbose):
        sys.stderr.write(msg)

def syscmd(cmd):
    log("running:\n\t"+cmd+"\n")
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % exitstatus
        sys.exit(10)
    log("... done\n")

def main(args):
    ## Parse arguments
    try:
        opts, args = getopt.getopt(args, "c:q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 0:
        print "incorrect number of arguments"
        usage(1)
    
    stream = sys.stdin
    
    countf = None
    global verbose
    for o, a in opts:
        if o == "-c":
            countf = a
        if o == "-q":
            verbose = False
    
    ## Compute stats for each network
    nodeCounts = dict()
    totNodes_list = []
    totLinks_list = []
    largestNodes_list = []
    largestLinks_list = []
    if verbose:
        print "id\ttotNodes\ttotLinks\tlargest_netNodes\tlargest_netLinks"
    for line in stream:
        line = line.rstrip("\n\r")
        if line.isspace():
            continue
        (forNodes, forInteractions) = mPathway.rSIF(line, "node")
        (revNodes, revInteractions) = mPathway.rSIF(line, "node", reverse = True)
        for i in forNodes.keys():
            if i not in nodeCounts:
                nodeCounts[i] = 0
            nodeCounts[i] += 1
        totalNodes = len(forNodes.keys())
        totalLinks = 0
        for i in forInteractions.keys():
            totalLinks += len(forInteractions[i].keys())
        (lNodes, lInteractions) = mPathway.largestConnected(forNodes, forInteractions, revInteractions)
        largestNodes = len(lNodes.keys())
        largestLinks = 0
        for i in lInteractions.keys():
            largestLinks += len(lInteractions[i].keys())
        totNodes_list.append(totalNodes)
        totLinks_list.append(totalLinks)
        largestNodes_list.append(largestNodes)
        largestLinks_list.append(largestLinks)
        print "%s\t%s\t%s\t%s\t%s" % (line, totalNodes, totalLinks, largestNodes, largestLinks)
    if (len(totNodes_list) > 1) & outStats:
        (mean_totNodes, std_totNodes) = mCalculate.mean_std(totNodes_list)
        (mean_totLinks, std_totLinks) = mCalculate.mean_std(totLinks_list)
        (mean_largestNodes, std_largestNodes) = mCalculate.mean_std(largestNodes_list)
        (mean_largestLinks, std_largestLinks) = mCalculate.mean_std(largestLinks_list)
        print "MEAN\t%s\t%s\t%s\t%s" % (mean_totNodes, mean_totLinks, mean_largestNodes, mean_largestLinks)
        print "STD\t%s\t%s\t%s\t%s" % (std_totNodes, std_totLinks, std_largestNodes, std_largestLinks)
    if countf != None:
        features = nodeCounts.keys()
        features.sort(lambda x, y: cmp(nodeCounts[y],nodeCounts[x]))
        f = open(countf, "w")
        for i in features:
            f.write("%s\t%s\n" % (i, nodeCounts[i]/float(len(totNodes_list))))
        f.close()

if __name__ == "__main__":
    main(sys.argv[1:])
