#!/usr/bin/env python
"""connect.py: 

Usage:
  connect.py [options]

Options:
  -q            run quietly
"""
## Written By: Sam Ng
## Last Updated: ##/##/####
import getopt, os, re, sys
import mData, mPathway

verbose = True
sessionName = "metabolism"
maxDistance = 9

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
        opts, args = getopt.getopt(args, "q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 3:
        log("ERROR: incorrect number of arguments", die = True)
    
    featureFile = args[0]
    pathwayFile = args[1]
    scoreFile = args[2]
    
    global verbose
    for o, a in opts:
        if o == "-q":
            verbose = False
    
    ## execute
    featureList = mData.rList(featureFile)
    (gNodes, gInteractions) = mPathway.rPathway(pathwayFile)
    scoreMap = {}
    scoreMap[sessionName] = mData.r2Col(scoreFile)
    
    ## find connected
    connectList = set()
    for source in featureList:
        if source not in gNodes:
            continue
        for target in featureList:
            if target not in gNodes:
                continue
            if source == target:
                continue
            paths = mPathway.shortestPath(source, target, gInteractions, maxDistance = maxDistance)
            if len(paths) == 0:
                log("%s /> %s\n" % (source, target))
            else:
                log("%s -> %s\n" % (source, target))
            for path in paths:
                for feature in path:
                    connectList.update([feature])
    connectList = list(connectList)
    (cNodes, cInteractions) = mPathway.constructInteractions(connectList, gNodes, gInteractions)
    
    ## output
    os.system("mkdir LAYOUT")
    os.system("mkdir LAYOUT/%s" % (sessionName))
    mPathway.wSIF("LAYOUT/%s/%s.sif" % (sessionName, sessionName), cInteractions)
    #mPathway.wNodeAttributes(gNodes, scoreMap = scoreMap, directory = "LAYOUT")
if __name__ == "__main__":
    main(sys.argv[1:])
