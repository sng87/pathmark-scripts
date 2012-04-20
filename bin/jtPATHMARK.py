#!/usr/bin/env python
"""jtPATHMARK.py:

Usage:
  jtPATHMARK.py

Optional Configuration Files:
  
"""
## Written By: Sam Ng
## Last Updated: 04/19/2012
import math, os, os.path, sys, random, re
from optparse import OptionParser
from copy import deepcopy
from mPATHMARK import *

from jobTree.src.bioio import logger
from jobTree.src.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

## executables and directories
lmExec = "lm.R"
occamExec = "OCCAM.py"
pathmarkExec = "PATHMARK.py"
statisticsExec = "statistics-subnets.py"

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
    def __init__(self, dataFile, sampleList, directory):
        Target.__init__(self, time=1000)
        self.dataFile = dataFile
        self.sampleList = sampleList
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## write data
        if not os.path.exists("real.tab"):
            rwCRSData("real.tab", "%s" % (self.dataFile), useCols = self.sampleList)

class jtNData(Target):
    def __init__(self, null, dataFile, sampleList, directory):
        Target.__init__(self, time=1000)
        self.null = null
        self.dataFile = dataFile
        self.sampleList = sampleList
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## permute features
        colMap = {}
        for sample in self.sampleList:
            rindex = random.randint(1,5)
            colMap["na_iter_%s_%s" % (rindex, sample)] = sample
            colMap[sample] = "na_iter_%s_%s" % (rindex, sample)
        
        ## write data
        if not os.path.exists("null_%s.tab" % (self.null)):
            rwCRSData("null_%s.tab" % (self.null), "%s" % (self.dataFile), colMap = colMap, useCols = self.sampleList)

class jtCmd(Target):
    def __init__(self, cmd, directory):
        Target.__init__(self, time=1000)
        self.cmd = cmd
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        system(self.cmd)

class prepareOCCAM(Target):
    def __init__(self, paradigmPathway, scoreFile, phenotypeFile, subtractFile, dataFile, sampleList, filterParams, nNulls, directory):
        Target.__init__(self, time=10000)
        self.paradigmPathway = paradigmPathway
        self.scoreFile = scoreFile
        self.phenotypeFile = phenotypeFile
        self.subtractFile = subtractFile
        self.dataFile = dataFile
        self.sampleList = sampleList
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## create real and null matricies from merge_merged.all.tab
        if self.dataFile is not None:
            self.addChildTarget(jtData(self.dataFile, self.sampleList, self.directory))
            for null in range(1, self.nNulls + 1):
                self.addChildTarget(jtNData(null, self.dataFile, self.sampleList, self.directory))
        self.setFollowOnTarget(runOCCAM(self.paradigmPathway, self.scoreFile, self.phenotypeFile, self.subtractFile, self.dataFile, self.sampleList, self.filterParams, self.nNulls, self.directory))
        
class runOCCAM(Target):
    def __init__(self, paradigmPathway, scoreFile, phenotypeFile, subtractFile, dataFile, sampleList, filterParams, nNulls, directory):
        Target.__init__(self, time=10000)
        self.paradigmPathway = paradigmPathway
        self.scoreFile = scoreFile
        self.phenotypeFile = phenotypeFile
        self.subtractFile = subtractFile
        self.dataFile = dataFile
        self.sampleList = sampleList
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## run lm.R or OCCAM.py
        if self.phenotypeFile is not None:
            phenotypeName = re.split("/", self.phenotypeFile)[-1].rstrip(".tab")
        else:
            phenotypeName = "feature"
        if self.subtractFile is not None:
            if not os.path.exists("OCCAM__%s__real" % (phenotypeName)):
                self.addChildTarget(jtCmd("%s real.tab %s %s" % (lmExec, self.phenotypeFile, self.subtractFile), self.directory))
            for null in range(1, self.nNulls + 1):
                if not os.path.exists("OCCAM__%s__null_%s" % (phenotypeName, null)):
                    self.addChildTarget(jtCmd("%s null_%s.tab %s %s" % (lmExec, null, self.phenotypeFile, self.subtractFile), self.directory))
        elif self.phenotypeFile is not None:
            if not os.path.exists("OCCAM__%s__real" % (phenotypeName)):
                self.addChildTarget(jtCmd("%s %s real.tab" % (occamExec, self.phenotypeFile), self.directory))
            for null in range(1, self.nNulls + 1):
                if not os.path.exists("OCCAM__%s__null_%s" % (phenotypeName, null)):
                    self.addChildTarget(jtCmd("%s %s null_%s.tab" % (occamExec, self.phenotypeFile, null), self.directory))
        self.setFollowOnTarget(branchPATHMARK(self.paradigmPathway, self.scoreFile, self.phenotypeFile, self.dataFile, self.sampleList, self.filterParams, self.nNulls, self.directory))

class branchPATHMARK(Target):
    def __init__(self, paradigmPathway, scoreFile, phenotypeFile, dataFile, sampleList, filterParams, nNulls, directory):
        Target.__init__(self, time=10000)
        self.paradigmPathway = paradigmPathway
        self.scoreFile = scoreFile
        self.phenotypeFile = phenotypeFile
        self.dataFile = dataFile
        self.sampleList = sampleList
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.directory = directory
    def run(self):
        layoutDir = "%s/LAYOUT" % (self.directory)
        if not os.path.exists(layoutDir):
            system("mkdir %s" % (layoutDir))
        os.chdir(layoutDir)
        
        ## link pathway
        if not os.path.exists(re.split("/", self.paradigmPathway)[-1]):
            if self.paradigmPathway.startswith("/"):
                system("ln -s %s ." % (self.paradigmPathway))
            else:
                system("ln -s ../%s ." % (self.paradigmPathway))
        
        ## link real_results.all.tab
        if self.phenotypeFile is not None:
            phenotypeName = re.split("/", self.phenotypeFile)[-1].rstrip(".tab")
        else:
            phenotypeName = "feature"
        if not os.path.exists("real_results.all.tab"):
            if self.scoreFile is None:
                system("ln ../OCCAM__%s__real/results.tab real_results.all.tab" % (phenotypeName))
            else:
                if self.scoreFile.startswith("/"):
                    system("ln %s real_results.all.tab" % (self.scoreFile))
                else:
                    system("ln ../%s real_results.all.tab" % (self.scoreFile))
        if len(retColumns("real_results.all.tab")) == 0:
            log("Differential analysis failed! Exiting PATHMARK ...\n", die = True)
        
        ## iterate through occamPhenotypes
        occamPhenotypes = retColumns("real_results.all.tab")
        for occamPhenotype in occamPhenotypes:
            self.addChildTarget(runPATHMARK(occamPhenotype, self.paradigmPathway, self.scoreFile, self.phenotypeFile, self.dataFile, self.sampleList, self.filterParams, self.nNulls, self.directory))
        self.setFollowOnTarget(cleanup(self.directory))

class runPATHMARK(Target):
    def __init__(self, occamPhenotype, paradigmPathway, scoreFile, phenotypeFile, dataFile, sampleList, filterParams, nNulls, directory):
        Target.__init__(self, time=10000)
        self.occamPhenotype = occamPhenotype
        self.paradigmPathway = paradigmPathway
        self.scoreFile = scoreFile
        self.phenotypeFile = phenotypeFile
        self.dataFile = dataFile
        self.sampleList = sampleList
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.directory = directory
    def run(self):
        layoutDir = "%s/LAYOUT" % (self.directory)
        os.chdir(layoutDir)
        
        ## aggregate null scores
        if self.phenotypeFile is not None:
            phenotypeName = re.split("/", self.phenotypeFile)[-1].rstrip(".tab")
        else:
            phenotypeName = "feature"
        if self.nNulls > 0:
            if not os.path.exists("null_results.%s.tab" % (self.occamPhenotype)):
                nullScores = {}
                for null in range(1, self.nNulls + 1):
                    if len(retColumns("../OCCAM__%s__null_%s/results.tab" % (phenotypeName, null))) == 0:
                        ## this is an error right now
                        continue
                    nullScores["N%s" % (null)] = rCRSData("../OCCAM__%s__null_%s/results.tab" % (phenotypeName, null))[self.occamPhenotype]
                wCRSData("null_results.%s.tab" % (self.occamPhenotype), nullScores)
        
        ## run pathmark
        system("%s -b %s -f %s -n real_results.all.tab >& %s.params" % (pathmarkExec, self.filterParams, self.occamPhenotype, self.occamPhenotype))
        if self.nNulls > 0:
            system("%s -b %s -s %s.params -d NULL_%s null_results.%s.tab" % (pathmarkExec, self.filterParams, self.occamPhenotype, self.occamPhenotype, self.occamPhenotype))
            self.setFollowOnTarget(backgroundPATHMARK(self.occamPhenotype, self.nNulls, self.directory))

class backgroundPATHMARK(Target):
    def __init__(self, occamPhenotype, nNulls, directory):
        Target.__init__(self, time=10000)
        self.occamPhenotype = occamPhenotype
        self.nNulls = nNulls
        self.directory = directory
    def run(self):
        layoutDir = "%s/LAYOUT" % (self.directory)
        os.chdir(layoutDir)
        
        system("ls %s/*_nodrug.sif | %s > stats_%s.tab" % (self.occamPhenotype, statisticsExec, self.occamPhenotype))
        system("ls NULL_%s/*_nodrug.sif | %s -c counts_NULL_%s.tab > stats_NULL_%s.tab" % (self.occamPhenotype, statisticsExec, self.occamPhenotype, self.occamPhenotype))
        system("../background.R %s" % (self.occamPhenotype))

class cleanup(Target):
    def __init__(self, directory):
        Target.__init__(self, time=10000)
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        system("rm -rf real* null* OCCAM__* background.R LAYOUT/*.params LAYOUT/real_results.* LAYOUT/null_results.*")

def main():
    ## parse arguments
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help="Add as a child of jobFile rather " +
                      "than making a new jobTree")
    parser.add_option("-b", "--boundaries", dest="boundParams", default="\"0.0;0.0\"")
    parser.add_option("-r", "--randoms", dest="nBackground", default="0")
    options, args = parser.parse_args()
    print "Using Batch System '" + options.batchSystem + "'"
    
    ## clean
    if len(args) == 1:
        if args[0] == "clean":
            print "rm -rf real* null* OCCAM__* LAYOUT background.R .jobTree"
            os.system("rm -rf real* null* OCCAM__* LAYOUT background.R .jobTree")
            sys.exit(0)
    
    ## parse arguments
    assert ((len(args) == 2) or (len(args) == 3))
    if len(args) == 2:
        paradigmPathway = args[0] 
        scoreFile = args[1]
        phenotypeFile = None
        dataFile = None
        sampleList = None
        filterParams = options.boundParams
        nNulls = 0
        assert(os.path.exists(paradigmPathway))
        assert(os.path.exists(scoreFile))
    else:
        paradigmPathway = args[0]
        scoreFile = None
        phenotypeFile = args[2]
        dataFile = args[1]
        sampleList = []
        for sample in retColumns(dataFile):
            if not sample.startswith("na_iter"):
                sampleList.append(sample)
        filterParams = options.boundParams
        nNulls = int(options.nBackground)
        assert(os.path.exists(paradigmPathway))
        assert(os.path.exists(phenotypeFile))
        assert(os.path.exists(dataFile))
    
    ## run
    logger.info("options: " + str(options))
    logger.info("starting make")
    writeScripts()
    s = Stack(prepareOCCAM(paradigmPathway, scoreFile, phenotypeFile, None, dataFile, sampleList, filterParams, nNulls, os.getcwd()))
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
