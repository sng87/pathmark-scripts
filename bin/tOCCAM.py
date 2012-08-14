#!/usr/bin/env python
"""tOCCAM.py: 

Usage:
  tOCCAM.py [options] phenotype data

Options:
  -q            run quietly
"""
## Written By: Sam Ng
## Last Updated: ##/##/####
import getopt, math, os, re, sys
import mData

verbose = True
phenIndex = 0

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

def mean_std(inList, sample = True):
    """Calculates mean and std"""
    cList = mData.floatList(inList)
    if len(cList) == 0:
        mean = "NA"
        std = "NA"
    else:
        mean = sum(cList)/float(len(cList))
        std = 0.0
        for i in cList:
            std += (i-mean)**2
        if len(cList) > 1:
            if sample:
                std = math.sqrt(std/(len(cList)-1))
            else:
                std = math.sqrt(std/len(cList))
        else:
            std = 0.0
    return(mean, std)

def ttest(values0, values1, alpha = 0.05):
    (mean0, std0) = mean_std(values0)
    (mean1, std1) = mean_std(values1)
    try:
        tval = (mean1-mean0)/(math.sqrt((std1**2)/len(values1)+(std0**2)/len(values0))+alpha)
    except ZeroDivisionError:
        tval = "NA"
    except TypeError:
        tval = "NA"
    return(tval)

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 2:
        log("ERROR: incorrect number of arguments", die = True)
    
    phenotypeFile = args[0]
    dataFile = args[1]
    
    global verbose
    for o, a in opts:
        if o == "-q":
            verbose = False
    
    ## execute
    phenotypeName = re.split("/", phenotypeFile)[-1].rstrip(".tab")
    dataName = re.split("/", dataFile)[-1].rstrip(".tab")
    outputDir = "OCCAM__%s__%s" % (phenotypeName, dataName)
    syscmd("mkdir %s" % (outputDir))
    (phenData, phenColumns, phenRows) = mData.rCRSData(phenotypeFile, retFeatures = True)
    (matData, matColumns, matRows) = mData.rCRSData(dataFile, retFeatures = True)
    
    ## samples
    posSamples = []
    negSamples = []
    for sample in phenRows:
        if sample not in matColumns:
            continue
        if phenData[phenColumns[phenIndex]][sample] == "+":
            posSamples.append(sample)
        elif phenData[phenColumns[phenIndex]][sample] == "-":
            negSamples.append(sample)
    
    ## output
    f = open("%s/results.tab" % (outputDir), "w")
    f.write("id\t%s\n" % (phenColumns[phenIndex]))
    for row in matRows:
        val = ttest([matData[i][row] for i in negSamples], [matData[i][row] for i in posSamples])
        f.write("%s\t%s\n" % (row, val))
    f.close()

if __name__ == "__main__":
    main(sys.argv[1:])
