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

basedir = os.path.dirname(os.path.abspath(__file__))

basepathway = os.path.join(basedir, "p_global_five3_v2.zip")

lmExec = os.path.join(basedir, "lm.R")
occamExec = os.path.join(basedir, "OCCAM.py")
ttestExec = os.path.join(basedir, "tOCCAM.py")
pathmarkExec = os.path.join(basedir, "PATHMARK.py")
statisticsExec = os.path.join(basedir, "statistics-subnets.py")

def writeScripts():
    """creates the R scripts necessary for plotting"""
    backgroundR = """#!/usr/bin/env Rscript
    args = commandArgs(TRUE)
    phenotype = args[1]
    Real = read.table(paste("stats_", phenotype, ".tab", sep=""), header=TRUE)
    Nulls = read.table(paste("stats_NULL_", phenotype, ".tab", sep=""), header=TRUE)
    nbreaks = 60
    
    zscore = c(as.character((Real$totNodes-mean(Nulls$totNodes))/sd(Nulls$totNodes)), as.character((Real$totLinks-mean(Nulls$totLinks))/sd(Nulls$totLinks)), as.character((Real$largest_netNodes-mean(Nulls$largest_netNodes))/sd(Nulls$largest_netNodes)), as.character((Real$largest_netLinks-mean(Nulls$largest_netLinks))/sd(Nulls$largest_netLinks)))
    fileConn = file(paste(phenotype, ".stats", sep=""))
    writeLines(zscore, fileConn)
    close(fileConn)
    
    xrange = c(min(Nulls$totNodes, Real$totNodes)-50, max(Nulls$totNodes, Real$totNodes)+50)
    png(paste(phenotype, "_total_netNodes.png", sep=""), heigh=720, width=1280)
    hist(Nulls$totNodes, breaks=nbreaks, xlim=xrange, xlab="Number", main=paste("Number of Nodes for Subnet, z = ", zscore[1], sep=""))
    abline(v = Real$totNodes, col="red", lty = 2)
    dev.off()
    
    xrange = c(min(Nulls$totLinks, Real$totLinks)-50, max(Nulls$totLinks, Real$totLinks)+50)
    png(paste(phenotype, "_total_netLinks.png", sep=""), heigh=720, width=1280)
    hist(Nulls$totLinks, breaks=nbreaks, xlim=xrange, xlab="Number", main=paste("Number of Links for Subnet, z = ", zscore[2], sep=""))
    abline(v = Real$totLinks, col="red", lty = 2)
    dev.off()
    
    xrange = c(min(Nulls$largest_netNodes, Real$largest_netNodes)-50, max(Nulls$largest_netNodes, Real$largest_netNodes)+50)
    png(paste(phenotype, "_largest_netNodes.png", sep=""), heigh=720, width=1280)
    hist(Nulls$largest_netNodes, breaks=nbreaks, xlim=xrange, xlab="Number", main=paste("Number of Nodes for Largest Component, z = ", zscore[3], sep=""))
    abline(v = Real$largest_netNodes, col="red", lty = 2)
    dev.off()
    
    xrange = c(min(Nulls$largest_netLinks, Real$largest_netLinks)-50, max(Nulls$largest_netLinks, Real$largest_netLinks)+50)
    png(paste(phenotype, "_largest_netLinks.png", sep=""), heigh=720, width=1280)
    hist(Nulls$largest_netLinks, breaks=nbreaks, xlim=xrange, xlab="Number", main=paste("Number of Links for Largest Component, z = ", zscore[4], sep=""))
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
    def __init__(self, paradigmPathway, scoreFile, phenotypeFile, subtractFile, dataFile, sampleList, filterParams, nNulls, outputZip, directory):
        Target.__init__(self, time=10000)
        self.paradigmPathway = paradigmPathway
        self.scoreFile = scoreFile
        self.phenotypeFile = phenotypeFile
        self.subtractFile = subtractFile
        self.dataFile = dataFile
        self.sampleList = sampleList
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.outputZip = outputZip
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## create real and null matricies from merge_merged.all.tab
        if self.dataFile is not None:
            self.addChildTarget(jtData(self.dataFile, self.sampleList, self.directory))
            for null in range(1, self.nNulls + 1):
                self.addChildTarget(jtNData(null, self.dataFile, self.sampleList, self.directory))
        self.setFollowOnTarget(runOCCAM(self.paradigmPathway, self.scoreFile, self.phenotypeFile, self.subtractFile, self.dataFile, self.sampleList, self.filterParams, self.nNulls, self.outputZip, self.directory))
        
class runOCCAM(Target):
    def __init__(self, paradigmPathway, scoreFile, phenotypeFile, subtractFile, dataFile, sampleList, filterParams, nNulls, outputZip, directory):
        Target.__init__(self, time=10000)
        self.paradigmPathway = paradigmPathway
        self.scoreFile = scoreFile
        self.phenotypeFile = phenotypeFile
        self.subtractFile = subtractFile
        self.dataFile = dataFile
        self.sampleList = sampleList
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.outputZip = outputZip
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## run lm.R or OCCAM.py
        if self.phenotypeFile is not None:
            phenotypeName = re.split("/", self.phenotypeFile)[-1]
        if self.subtractFile is not None:
            if not os.path.exists("OCCAM__%s__real.tab" % (phenotypeName)):
                self.addChildTarget(jtCmd("%s real.tab %s %s" % (lmExec, self.phenotypeFile, self.subtractFile), self.directory))
            for null in range(1, self.nNulls + 1):
                if not os.path.exists("OCCAM__%s__null_%s.tab" % (phenotypeName, null)):
                    self.addChildTarget(jtCmd("%s null_%s.tab %s %s" % (lmExec, null, self.phenotypeFile, self.subtractFile), self.directory))
        elif self.phenotypeFile is not None:
            if not os.path.exists("OCCAM__%s__real.tab" % (phenotypeName)):
                self.addChildTarget(jtCmd("%s %s %s real.tab" % (sys.executable, occamExec, self.phenotypeFile), self.directory))
            for null in range(1, self.nNulls + 1):
                if not os.path.exists("OCCAM__%s__null_%s.tab" % (phenotypeName, null)):
                    self.addChildTarget(jtCmd("%s %s %s null_%s.tab" % (sys.executable, occamExec, self.phenotypeFile, null), self.directory))
        self.setFollowOnTarget(branchPATHMARK(self.paradigmPathway, self.scoreFile, self.phenotypeFile, self.dataFile, self.sampleList, self.filterParams, self.nNulls, self.outputZip, self.directory))

class branchPATHMARK(Target):
    def __init__(self, paradigmPathway, scoreFile, phenotypeFile, dataFile, sampleList, filterParams, nNulls, outputZip, directory):
        Target.__init__(self, time=10000)
        self.paradigmPathway = paradigmPathway
        self.scoreFile = scoreFile
        self.phenotypeFile = phenotypeFile
        self.dataFile = dataFile
        self.sampleList = sampleList
        self.filterParams = filterParams
        self.nNulls = nNulls
        self.outputZip = outputZip
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
        if not os.path.exists("real_results.all.tab"):
            if self.scoreFile is None:
                phenotypeName = re.split("/", self.phenotypeFile)[-1]
                system("ln ../OCCAM__%s__real.tab/results.tab real_results.all.tab" % (phenotypeName))
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
        self.setFollowOnTarget(cleanup(self.outputZip, self.directory))

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
        if self.nNulls > 0:
            phenotypeName = re.split("/", self.phenotypeFile)[-1]
            if not os.path.exists("null_results.%s.tab" % (self.occamPhenotype)):
                nullScores = {}
                for null in range(1, self.nNulls + 1):
                    if len(retColumns("../OCCAM__%s__null_%s.tab/results.tab" % (phenotypeName, null))) == 0:
                        ## this is an error right now
                        continue
                    nullScores["N%s" % (null)] = rCRSData("../OCCAM__%s__null_%s.tab/results.tab" % (phenotypeName, null))[self.occamPhenotype]
                wCRSData("null_results.%s.tab" % (self.occamPhenotype), nullScores)
        
        ## run pathmark
        system("%s %s -l %s.params -b \"%s\" -f %s -n real_results.all.tab" % (sys.executable, pathmarkExec, self.occamPhenotype, self.filterParams, self.occamPhenotype))
        if self.nNulls > 0:
            system("%s %s -b \"%s\" -s %s.params -d NULL_%s null_results.%s.tab" % (sys.executable, pathmarkExec, self.filterParams, self.occamPhenotype, self.occamPhenotype, self.occamPhenotype))
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
        
        system("ls %s/*_nodrug.sif | %s %s > stats_%s.tab" % (self.occamPhenotype, sys.executable, statisticsExec, self.occamPhenotype))
        system("ls NULL_%s/*_nodrug.sif | %s %s -c counts_NULL_%s.tab > stats_NULL_%s.tab" % (self.occamPhenotype, sys.executable, statisticsExec, self.occamPhenotype, self.occamPhenotype))
        system("../background.R %s" % (self.occamPhenotype))

class cleanup(Target):
    def __init__(self, outputZip, directory):
        Target.__init__(self, time=10000)
        self.outputZip = outputZip
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        system("rm -rf real* null* OCCAM__* background.R LAYOUT/*.params LAYOUT/real_results.* LAYOUT/null_results.* LAYOUT/*.tab LAYOUT/NULL_*")
        if self.outputZip is not None:
            system("zip -r LAYOUT.zip LAYOUT")
            system("mv -f LAYOUT.zip %s" % (self.outputZip))

def main():
    ## parse arguments
    parser = OptionParser(usage = "%prog [options] network IPL-matrix features")
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help="Add as a child of jobFile rather " +
                      "than making a new jobTree")
    parser.add_option("-w", "--workdir", dest="workdir", help="Common Work directory", default="./")
    parser.add_option("-i", "--ipl", dest="iplFile", default = None)
    parser.add_option("-p", "--pathway", dest="pathwayZip", default=None)
    parser.add_option("-c", "--phenotype", dest="phenotypeFile", default=None)
    parser.add_option("-o", "--oz", dest="outputZip", default=None)
    parser.add_option("-s", "--score", dest="scoreFile", default=None)
    parser.add_option("-f", "--filter", dest="filterParams", default="0.0;0.0")
    parser.add_option("-b", "--background", dest="nBackground", default="0")
    options, args = parser.parse_args()
    print "Using Batch System '" + options.batchSystem + "'"
    
    ## clean
    if len(args) == 1:
        if args[0] == "clean":
            print "rm -rf real* null* OCCAM__* LAYOUT background.R .jobTree"
            system("rm -rf real* null* OCCAM__* LAYOUT background.R .jobTree")
            sys.exit(0)
    
    ## parse arguments
    assert ((len(args) == 0) or (len(args) == 2) or (len(args) == 3))
    if len(args) == 0:
        pathwayZip = options.pathwayZip if options.pathwayZip is not None else basepathway
        pathwayLib = os.path.join(options.workdir, "pathway")
        system("unzip %s -d %s" % (pathwayZip, pathwayLib))
        paradigmPathway = None
        for file in os.listdir(pathwayLib):
            if file.endswith("_pathway.tab"):
                paradigmPathway = "%s/%s" % (pathwayLib, file)
                break
        scoreFile = None
        phenotypeFile = options.phenotypeFile
        dataFile = options.iplFile
        sampleList = []
        for sample in retColumns(dataFile):
            if not sample.startswith("na_iter"):
                sampleList.append(sample)
        filterParams = options.filterParams
        nNulls = int(options.nBackground)
        outputZip = options.outputZip
        assert(os.path.exists(paradigmPathway))
        assert(os.path.exists(phenotypeFile))
        assert(os.path.exists(dataFile))
    elif len(args) == 2:
        paradigmPathway = args[0] 
        scoreFile = args[1]
        phenotypeFile = None
        dataFile = None
        sampleList = None
        filterParams = options.filterParams
        nNulls = 0
        outputZip = options.outputZip
        assert(os.path.exists(paradigmPathway))
        assert(os.path.exists(scoreFile))
    elif len(args) == 3:
        paradigmPathway = args[0]
        scoreFile = None
        phenotypeFile = args[2]
        dataFile = args[1]
        sampleList = []
        for sample in retColumns(dataFile):
            if not sample.startswith("na_iter"):
                sampleList.append(sample)
        filterParams = options.filterParams
        nNulls = int(options.nBackground)
        outputZip = options.outputZip
        assert(os.path.exists(paradigmPathway))
        assert(os.path.exists(phenotypeFile))
        assert(os.path.exists(dataFile))
    
    ## run
    logger.info("options: " + str(options))
    logger.info("starting make")
    writeScripts()
    s = Stack(prepareOCCAM(paradigmPathway, scoreFile, phenotypeFile, None, dataFile, sampleList, filterParams, nNulls, outputZip, os.getcwd()))
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
