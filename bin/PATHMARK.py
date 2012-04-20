#!/usr/bin/env python
"""PATHMARK.py: identifies subnets in a score-matrix (feature x phenotype)

Usage:
  PATHMARK.py [options] scoref [scoref ...]

Options:
  -p  str                     pathway file to superimpose scores onto (default: ./*_pathway.tab)
  -f  str[,str,...]           features to run PATHMARK on (default: all)
  -s  flt;flt[,flt;flt,...]   force mean and standard deviation statistics for score files
  -b  flt;flt                 filter lower and upper boundary parameters (default: 0;0)
  -d  str                     output to user-specified directory (default: feature/)
  -n                          output node attributes for cytoscape
  -t                          output paradigm net
  -q                          run quietly
"""
## Written By: Sam Ng
## Last Updated: 7/23/11
import os, os.path, sys, getopt, re
from mPATHMARK import *
from copy import deepcopy

verbose = True
outputCleaned = False

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg, die = False):
    if (verbose):
        sys.stderr.write(msg)
    if (die):
        sys.exit(1)

def syscmd(cmd):
    log("running: "+cmd)
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % exitstatus
        sys.exit(10)
    log(" ... done\n")

def PATHMARK(files, globalPathway, features = None, statLine = None, 
             filterBounds = [0,0], outDir = None, outputAttributes = False, 
             outputPARADIGM = False, selectionRule = "OR", topDisconnected = 100):
    filterString = "%s_%s" % (filterBounds[0], filterBounds[1])
    
    ## read global pathway
    (gNodes, gInteractions) = rPathway(globalPathway)
    
    ## read scores
    uData = {}
    sData = {}
    for i in range(len(files)):
        uData[i] = rCRSData(files[i])
        sData[i] = {}
        for j in uData[i].keys():
            sData[i][j] = {}
            for k in uData[i][j].keys():
                try:
                    sData[i][j][k] = abs(float(uData[i][j][k]))
                except ValueError:
                    sData[i][j][k] = "NA"
    
    ## iterate features
    for feature in sData[0].keys():
        if features is not None:
            if feature not in features:
                continue
        pNodes = {}
        pInteractions = {}
        
        ## compute score statistics
        pStats = []
        if statLine is None:
            for i in range(len(sData.keys())):
                pStats.append(mean_std(sData[i][feature].values()))
        else:
            for i in re.split(",",statLine):
                (v1, v2) = re.split(";",i)
                pStats.append((float(v1), float(v2)))
        
        ## log statistics
        log("%s\t%s;%s" % (feature, pStats[0][0], pStats[0][1]))
        for i in range(1, len(pStats)):
            log(",%s;%s" % (pStats[i][0], pStats[i][1]))
        log("\n")
        
        ## add top scoring links
        for source in gInteractions.keys():
            if source not in sData[0][feature]:
                continue
            elif sData[0][feature][source] == "NA":
                continue
            for target in gInteractions[source].keys():
                if target not in sData[0][feature]:
                    continue
                elif sData[0][feature][target] == "NA":
                    continue
                
                if selectLink(feature, source, target, sData, pStats, filterBounds, selectionRule = selectionRule):
                    (pNodes, pInteractions) = addLink(source, target, pNodes, pInteractions, gNodes, gInteractions)
        
        ## add top scoring disconnected nodes
        sortedTop = []
        for node in sData[0][feature].keys():
            if node not in gNodes:
                continue
            if gNodes[node] in ["protein"]:
                sortedTop.append(node)
        sortedTop.sort(lambda x, y: cmp(sData[0][feature][y],sData[0][feature][x]))
        while (sData[0][feature][sortedTop[0]] == "NA"):
            sortedTop.pop(0)
            if len(sortedTop) == 0:
                break
        for i in range(topDisconnected):
            if i > len(sortedTop)-1:
                break
            if sData[0][feature][sortedTop[i]] < pStats[0][0]+filterBounds[0]*pStats[0][1]:
                break
            if sortedTop[i] not in pNodes:
                pNodes[sortedTop[i]] = gNodes[sortedTop[i]]
                pInteractions[sortedTop[i]] = {}
                pInteractions[sortedTop[i]]["__DISCONNECTED__"] = "-disconnected-"
        
        ## output node attributes
        if outputAttributes:
            wNodeAttributes(feature, gNodes, uData[0])
        
        ## output networks
        if outDir == None:
            wrtDir = feature
        else:
            wrtDir = outDir
        if not os.path.exists(wrtDir):
            syscmd("mkdir %s" % (wrtDir))
        (cpNodes, cpInteractions) = filterComplexesByGeneSupport(pNodes, pInteractions,
                                        reverseInteractions(pInteractions), gNodes,
                                        getComponentMap(gNodes, reverseInteractions(gInteractions)))
        if outputPARADIGM:
            wPARADIGM("%s/%s_%s_nodrug_pathway.tab" % (wrtDir, feature, filterString), pNodes, pInteractions)
            if outputCleaned:
                wPARADIGM("%s/%s_%s_nodrug_cleaned_pathway.tab" % (wrtDir, feature, filterString), cpNodes, cpInteractions)
        else:
            wSIF("%s/%s_%s_nodrug.sif" % (wrtDir, feature, filterString), pInteractions)
            if outputCleaned:
                wSIF("%s/%s_%s_nodrug_cleaned.sif" % (wrtDir, feature, filterString), cpInteractions)

if __name__ == "__main__":
    ## parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:f:s:b:d:ntq")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) < 1:
        print "incorrect number of arguments"
        usage(1)
    
    globalPathway = None
    features = None
    statLine = None
    filterBounds = [0,0]
    outDir = None
    outputAttributes = False
    outputPARADIGM = False
    selectionRule = "OR"
    for o, a in opts:
        if o == "-p":
            globalPathway = a
        elif o == "-f":
            features = re.split(";", a)
        elif o == "-s":
            statLine = a
            if os.path.exists(statLine):
                f = open(statLine, "r")
                statLine = re.split("\t", f.readline().rstrip("\r\n"))[1]
                f.close()
        elif o == "-b":
            (v1, v2) = re.split(";", a)
            filterBounds = [float(v1), float(v2)]
            filterBounds.sort()
        elif o == "-d":
            outDir = a
            if outDir.endswith("/"):
                outDir = outDir.rstrip("/")
        elif o == "-n":
            outputAttributes = True
        elif o == "-t":
            outputPARADIGM = True
        elif o == "-q":
            verbose = False
    if globalPathway is None:
        for file in os.listdir("."):
            if file.endswith("_pathway.tab"):
                globalPathway = file
    assert globalPathway is not None
    
    ## run
    PATHMARK(args, globalPathway, features = features, statLine = statLine, 
             filterBounds = filterBounds, outDir = outDir, outputAttributes = outputAttributes, 
             outputPARADIGM = outputPARADIGM, selectionRule = selectionRule, topDisconnected = 100)
