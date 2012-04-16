#!/usr/bin/env python
"""jtPATHMARK.py:

Usage:
  jtPATHMARK.py

Optional Configuration Files:
  
"""
## Written By: Sam Ng
## Last Updated: ##/##/####
import math, os, os.path, sys, random, re
from optparse import OptionParser
from copy import deepcopy
from mPATHMARK import *

from jobTree.src.bioio import logger
from jobTree.src.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

## default variables
paradigmPathway = "/hive/users/sng/map/pathwayFiles/global_five3_v2/pid_110725_pathway.tab"

## executables and directories
lmExec = "lm.R"
occamExec = "OCCAM.py"
pathmarkExec = "PATHMARK.py"
statisticsExec = "statistics-subnets.py"

## check files
assert os.path.exists("feature.tab")
assert os.path.exists("merge_merged.all.tab")
featureFile = "feature.tab"
occamPhenotype = retColumns("feature.tab")[0]
dataFile = "merge_merged.all.tab"
if os.path.exists("include.samples"):
    allSamples = rList("include.samples")
else:
    allSamples = []
    for sample in retColumns("merge_merged.all.tab"):
        if not sample.startswith("na_iter"):
            allSamples.append(sample)

def writeScripts():
    """creates the R scripts necessary for plotting"""
    backgroundR = """#!/usr/bin/env Rscript
    args = commandArgs(TRUE)
    phenotype = args[1]
    Real = read.table(paste("stats_", phenotype, ".tab", sep=""), header=TRUE)
    Nulls = read.table(paste("stats_NULL_", phenotype, ".tab", sep=""), header=TRUE)
    nbreaks = 60

    xrange = c(min(Nulls$totNodes, Real$totNodes)-50, max(Nulls$totNodes, Real$totNodes)+50)
    pdf(paste(phenotype, "_totNodes.pdf", sep=""), heigh=10, width=20)
    hist(Nulls$totNodes, breaks=nbreaks, xlim=xrange, xlab="Number", main=paste("Number of Nodes for Subnet, z = ", as.character((Real$totNodes-mean(Nulls$totNodes))/sd(Nulls$totNodes)), sep=""))
    abline(v = Real$totNodes, col="red", lty = 2)
    dev.off()

    xrange = c(min(Nulls$totLinks, Real$totLinks)-50, max(Nulls$totLinks, Real$totLinks)+50)
    pdf(paste(phenotype, "_totLinks.pdf", sep=""), heigh=10, width=20)
    hist(Nulls$totLinks, breaks=nbreaks, xlim=xrange, xlab="Number", main=paste("Number of Links for Subnet, z = ", as.character((Real$totLinks-mean(Nulls$totLinks))/sd(Nulls$totLinks)), sep=""))
    abline(v = Real$totLinks, col="red", lty = 2)
    dev.off()

    xrange = c(min(Nulls$largest_netNodes, Real$largest_netNodes)-50, max(Nulls$largest_netNodes, Real$largest_netNodes)+50)
    pdf(paste(phenotype, "_largest_netNodes.pdf", sep=""), heigh=10, width=20)
    hist(Nulls$largest_netNodes, breaks=nbreaks, xlim=xrange, xlab="Number", main=paste("Number of Nodes for Largest Component, z = ", as.character((Real$largest_netNodes-mean(Nulls$largest_netNodes))/sd(Nulls$largest_netNodes)), sep=""))
    abline(v = Real$largest_netNodes, col="red", lty = 2)
    dev.off()

    xrange = c(min(Nulls$largest_netLinks, Real$largest_netLinks)-50, max(Nulls$largest_netLinks, Real$largest_netLinks)+50)
    pdf(paste(phenotype, "_largest_netLinks.pdf", sep=""), heigh=10, width=20)
    hist(Nulls$largest_netLinks, breaks=nbreaks, xlim=xrange, xlab="Number", main=paste("Number of Links for Largest Component, z = ", as.character((Real$largest_netLinks-mean(Nulls$largest_netLinks))/sd(Nulls$largest_netLinks)), sep=""))
    abline(v = Real$largest_netLinks, col="red", lty = 2)
    dev.off()
    """
    
    f = open("background.R", "w")
    f.write(backgroundR)
    f.close
    system("chmod 755 *.R")

class jtData(Target):
    def __init__(self, directory):
        Target.__init__(self, time=1000)
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## write data
        if not os.path.exists("real.tab"):
            rwCRSData("real.tab", "%s" % (dataFile), useCols = allSamples)

class jtNData(Target):
    def __init__(self, null, directory):
        Target.__init__(self, time=1000)
        self.null = null
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## permute features
        colMap = {}
        for sample in allSamples:
            rindex = random.randint(1,5)
            colMap["na_iter_%s_%s" % (rindex, sample)] = sample
            colMap[sample] = "na_iter_%s_%s" % (rindex, sample)
        
        ## write data
        if not os.path.exists("null_%s.tab" % (self.null)):
            rwCRSData("null_%s.tab" % (self.null), "%s" % (dataFile), colMap = colMap, useCols = allSamples)

class jtCmd(Target):
    def __init__(self, cmd, directory):
        Target.__init__(self, time=1000)
        self.cmd = cmd
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        system(self.cmd)

class prepareOCCAM(Target):
    def __init__(self, filterParams, nNulls, directory):
        Target.__init__(self, time=10000)
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        self.addChildTarget(jtData(self.directory))
        for null in range(1, self.nNulls + 1):
            self.addChildTarget(jtNData(null, self.directory))
        self.setFollowOnTarget(runOCCAM(self.filterParams, self.nNulls, self.directory))
        
class runOCCAM(Target):
    def __init__(self, filterParams, nNulls, directory):
        Target.__init__(self, time=10000)
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        if os.path.exists("subtract.tab"):
            if not os.path.exists("OCCAM__feature__real"):
                self.addChildTarget(jtCmd("%s real.tab %s subtract.tab" % (lmExec, featureFile), self.directory))
            for null in range(1, self.nNulls + 1):
                if not os.path.exists("OCCAM__feature__null_%s" % (null)):
                    self.addChildTarget(jtCmd("%s null_%s.tab %s subtract.tab" % (lmExec, null, featureFile), self.directory))
        else:
            if not os.path.exists("OCCAM__feature__real"):
                self.addChildTarget(jtCmd("%s %s real.tab" % (occamExec, featureFile), self.directory))
            for null in range(1, self.nNulls + 1):
                if not os.path.exists("OCCAM__feature__null_%s" % (null)):
                    self.addChildTarget(jtCmd("%s %s null_%s.tab" % (occamExec, featureFile, null), self.directory))
        self.setFollowOnTarget(runPATHMARK(self.filterParams, self.nNulls, self.directory))

class runPATHMARK(Target):
    def __init__(self, filterParams, nNulls, directory):
        Target.__init__(self, time=10000)
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.directory = directory
    def run(self):
        layoutDir = "%s/LAYOUT" % (self.directory)
        if not os.path.exists(layoutDir):
            system("mkdir %s" % (layoutDir))
        os.chdir(layoutDir)
        
        if not os.path.exists(re.split("/", paradigmPathway)[-1]):
            system("ln -s %s ." % (paradigmPathway))
        if not os.path.exists("real_results.tab"):
            system("ln ../OCCAM__feature__real/results.tab real_results.tab")
        if self.nNulls > 0:
            if not os.path.exists("null_results.tab"):
                for null in range(1, self.nNulls + 1):
                    if not os.path.exists("null_results.tmp"):
                        system("cat ../OCCAM__feature__null_%s/results.tab | transpose.pl | head -n1 >> null_results.tmp" % (null))
                    system("cat ../OCCAM__feature__null_%s/results.tab | transpose.pl | grep \"%s\" | sed 's/%s/N%s/' >> null_results.tmp" % (null, occamPhenotype, occamPhenotype, null))
                system("cat null_results.tmp | transpose.pl > null_results.tab")
                system("rm -f null_results.tmp")
        
        system("%s -b %s -f %s -n real_results.tab >& a1" % (pathmarkExec, self.filterParams, occamPhenotype))
        if self.nNulls > 0:
            system("%s -b %s -s a1 -d NULL_%s null_results.tab" % (pathmarkExec, self.filterParams, occamPhenotype))
            self.setFollowOnTarget(backgroundPATHMARK(self.nNulls, self.directory))

class backgroundPATHMARK(Target):
    def __init__(self, nNulls, directory):
        Target.__init__(self, time=10000)
        self.nNulls = nNulls
        self.directory = directory
    def run(self):
        layoutDir = "%s/LAYOUT" % (self.directory)
        os.chdir(layoutDir)
        
        system("ls %s/*_nodrug.sif | %s > stats_%s.tab" % (occamPhenotype, statisticsExec, occamPhenotype))
        system("ls NULL_%s/*_nodrug.sif | %s -c counts_NULL_%s.tab > stats_NULL_%s.tab" % (occamPhenotype, statisticsExec, occamPhenotype, occamPhenotype))
        system("../background.R %s" % (occamPhenotype))

def main():
    ## parse arguments
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help="Add as a child of jobFile rather " +
                      "than making a new jobTree")
    parser.add_option("-b", "--boundaries", dest="boundParams", default="\"1.0;1.0\"")
    parser.add_option("-n", "--nulls", dest="nNulls", default="0")
    options, args = parser.parse_args()
    print "Using Batch System '" + options.batchSystem + "'"
    
    filterParams = options.boundParams
    nNulls = int(options.nNulls)
    
    ## clean
    if len(args) == 1:
        if args[0] == "clean":
            print "rm -rf null* real* LAYOUT OCCAM* .jobTree"
            os.system("rm -rf null* real* LAYOUT OCCAM* .jobTree")
            sys.exit(0)
    
    ## run
    assert len(args) == 0
    logger.info("options: " + str(options))
    logger.info("starting make")
    writeScripts()
    s = Stack(prepareOCCAM(filterParams, nNulls, os.getcwd()))
    if options.jobFile:
        s.addToJobFile(options.jobFile)
    else:
        if options.jobTree == None:
            options.jobTree = "./.jobTree"
        
        failed = s.startJobTree(options)
        if failed:
            print ("%d jobs failed" % failed)
        else:
            logger.info("Run complete!")
            system("rm -rf .lastjobTree")
            system("mv .jobTree .lastjobTree")

if __name__ == "__main__":
    from jtPATHMARK import *
    if os.path.exists(".jobTree"):
        print "WARNING: .jobTree directory already exists"
    main()
