#!/usr/bin/env python

#import pdb
import math
import os
import os.path
import time

import re
import sys
import operator
import random
import shutil

import subprocess
from operator import itemgetter


def die(arg):
   if arg:
       sys.stderr.write(arg +"\n")
   sys.stderr.write("""Usage: OCCAM phenotypeFile.tab featureFile.tab
   in phenotypeFile, phenotypes are stored in columns, sample names are stored in rows
   in featureFile.tab, feature names are stored in rows, sample names are stored in columns

Author
     Ted Goldstein (c) 2011, all rights reserved. Permission to use and modify is granted under a Creative Commons Attribution licenses.
     See http://creativecommons.org/licenses/by/3.0/
""")
   raise Exception(arg)


SHORTCIRCUIT = 999
ShortCircuit = SHORTCIRCUIT

ProcessMap = {}



TOOLSDIRS = [ os.path.dirname( os.path.abspath(__file__) ) ]
if "OCCAM_TOOLS" in os.environ:
    TOOLSDIR.insert(0, os.environ["OCCAM_TOOLS"])

TOOLS = None
for  x in TOOLSDIRS:
   if os.path.exists(x):
       TOOLS = x
       break

if TOOLS == None:
    die("Cannot find TOOLS directory; expected one of the following directories:\n" + "\n".join(TOOLSDIRS) +"\n you can also set the environment varialble OCCAM_TOOLS")



class DichotomizedPhenotypes:

    def __init__(this, filename):
        MINCOUNT = 1
        Data = []
        headers = None
        allRanges=[]
        this.allRanges = allRanges
        allColumns = []
        newColumns = []
        plusMinusSet = set(['+', '-'])


        def fix(s):
             s = s.replace('/', "+")
             s = s.replace('"', "")
             s = s.replace(" ", "_")
             s = s.replace("(", "_")
             s = s.replace(")", "_")
             s = s.replace("<", "_")
             s = s.replace(">", "_")
             s = s.replace("&", "_")
             return s


        for line in open(filename):
            if line[-2] == "\r":
                words = line[:-2].split("\t")
            else:
                words = line[:-1].split("\t")
            for i,w in enumerate(words):
               if len(w) >= 3 and ((w[0] == '"' and  w[-1] == '"') or (w[0] == "'" and  w[-1] == "'")):
                  words[i] = w[1:-1]

            if headers == None:
                headers=words
                for i in xrange(0,len(headers)): 
                     headers[i] = fix(headers[i])
                     allRanges.append(set())
                     allColumns.append(list())
            else:
                 for i,w in enumerate(words):
                      if w == '' or w == "null": 
                          w = ''
                      else:
                          num  = re.match("[+-]?\d*(\.\d+)?([Ee][+-]?\d+)?", w).group(0)
                          if len(num) > 0:
                               w = num
                      if w != '':
                          allRanges[i].add(w)
                      allColumns[i].append(w)


        NewColumns = []
        set_0_1 = set(["0","1"])



        for i,range in enumerate(allRanges):
           if i > 0 and len(range) > 1:
              try:
                 count = 0
                 total = 0
                 for value in allColumns[i]:
                    if len(value) > 0 and value != "NA":
                        count = count + 1
                        total = total + float(value)
                 average = total / count

                 n = [fix(headers[i]+"_MEAN")]
                 if allRanges[i] == set_0_1:
                     raise Exception('set 0 1')
                 neg= pos = 0
                 for j,value in enumerate(allColumns[i]):
                    if len(value) > 0 and value != "NA":
                        if float(value) >= average:
                            n.append("+")
                            pos = pos + 1
                        elif float(value) < average:
                            n.append("-")
                            neg = neg + 1
                        else:
                            #pdb.set_trace()
                            pass
                    else:
                        n.append("")
                 if pos >= MINCOUNT and neg >= MINCOUNT:
                     newColumns.append(n)
                 else:
                     sys.stdout.write(str(n[0])+ " insufficient count pos="+str(pos)+ " neg="+str(neg)+"\n")
              except Exception:
                 if len(range) < 2 or len(range) > 10:
                     pass
                 elif len(range) == 2:
                    element = None
                    for preferred in ["+", "YES", "Y", "1", "True", "TRUE"]:
                       if preferred in range:
                          element = preferred
                    if element == None:
                        element = list(range).pop()

                    if allRanges[i] == plusMinusSet:
                        n = [fix(headers[i])]
                    else:
                        if len(headers[i]) < 3:
                            print >>log, headers[i],element
                        n = [fix(headers[i] + "=" +element+ "_PHENOTYPE")]
                    neg= pos = 0
                    for value in allColumns[i]:
                        if value in range:
                            if value == element:
                                n.append("+")
                                pos = pos + 1
                            else:
                                n.append("-")
                                neg = neg + 1
                        else:
                           n.append("")
                    if pos >= MINCOUNT and neg >= MINCOUNT:
                        newColumns.append(n)
                    else:
                        sys.stderr.write(str(n[0])+ " insufficient count pos="+str(pos)+ " neg="+str(neg)+"\n")
                 else:
                    for element in range:
                        if allRanges[i] == plusMinusSet:
                            n = [fix(headers[i])]
                        else:
                            if len(headers[i]) < 3:
                                print >>log, headers[i],element
                            n = [fix(headers[i] + "=" +element+ "_PHENOTYPE")]
                        pos = neg = 0
                        for value in allColumns[i]:
                            if value in range:
                                if value == element:
                                    n.append("+")
                                    pos = pos + 1
                                else:
                                    n.append("-")
                                    neg = neg + 1
                            else:
                               n.append("")
                        if pos >= MINCOUNT and neg >= MINCOUNT:
                            newColumns.append(n)
                        else:
                            sys.stderr.write(str(n[0])+ " insufficient count pos="+str(pos)+ " neg="+str(neg)+"\n")


        q = allColumns[0]
        q.insert(0, "id")
        newColumns.insert(0, q)
        this.newColumns = newColumns


        this.rowLabelsArray = []
        this.rowMap = {}
        for i,r in enumerate(newColumns[0]):
            if i > 0:
                this.rowLabelsArray.append(r)
                this.rowMap[r] = i

        this.columnLabelsArray = []
        this.columnMap = {}
        for i,c in enumerate(newColumns):
            l = c[0]
            if i > 0:
                this.columnLabelsArray.append(l)
                this.columnMap[l] = i


    def dump(this, file):
        for i,name in enumerate(this.newColumns[0]):
            file.write(name)
            for column in this.newColumns[1:]:
                file.write("\t")
                file.write(column[i])
            file.write("\n")


    def phenotypeNames(this): # phenotypes are stored in columns
       return this.phenotype.columnLabels()

    def rowLabels(this): 
       return this.rowLabelsArray

    def columnLabels(this): 
       return this.columnLabelsArray

    def valueAtRowColumn(this, row, column):
        return this.newColumns[this.columnMap[column]][this.rowMap[row]]

class NullCracker:
    def __init__(this, featureTable):

        columnCount = len(featureTable.cellNames())
        groupCount = columnCount / 6
        # ensure that we have a multiple of six
        assert groupCount*6 == columnCount
        menu=[]

        # Make groupCount groups of 5
        for i in xrange(0, groupCount):
           menu.append([])

        columns = featureTable.columns()
        k = 1
        for j in xrange(0,5):
            for i in xrange(0, groupCount):
                c = columns[k]
                menu[i].append(c)
                k = k + 1

        labels = featureTable.columnLabels()
        this.columnMap = {}
        for i, five in enumerate(menu):
           this.columnMap[this.canon(labels[i])] = five

        this.concepts = featureTable.featureNames()

    # this.canonicalize the sample name back to the original 
    def canon(this, s):
       m = re.match("na_iter_([1-5])_(.*)", s)
       if m:
           s = m.group(2)
       return s

    def processPhenotype(this, phenotypeName, hypotheses, nullHypotheses):
        samples= {}
        samplesDichotomy = []
       
        print >>log, "Null processing for:", phenotypeName
        nullPhen = nullHypotheses+"/"+phenotypeName
        if not os.path.exists(nullPhen): os.mkdir(nullPhen)

        dich = []
        order = []
        for line in open(hypotheses+"/"+phenotypeName+"/sam.report"):
            if line[-2] == "\r":
                words = line[:-2].split("\t")
            else:
                words = line[:-1].split("\t")
            sampleName = this.canon(words[0])
            order.append(this.columnMap[sampleName])
            dich.append(words[2])

        crops = []
        for i in xrange(0,SHORTCIRCUIT):
            if ((i+1)%500) == 0: print phenotypeName,(i+1)
            dir = nullPhen +"/"+str(i) 
            print >> log, dir
            if not os.path.exists(dir): os.mkdir(dir)

            permutation = []
            for column in order:
                f = random.sample(column, 1)
                assert len(f) == 1
                permutation.append(f[0])

            f = dir+"/sam.input"
            samInput = open(f,"w")
            # write the header
            for s in dich:
               samInput.write("\t")
               samInput.write(s)
            samInput.write("\n")

            # write the data
            for j in range(1, len(this.concepts)):
                samInput.write(this.concepts[j])
                for p in permutation:
                    samInput.write("\t")
                    samInput.write(p[j])
                samInput.write("\n")
            samInput.close()
            runSam(str(i), dir)
         

def runSam(key, dir):
    print >>para, "para create",dir
    samCmd = dir+"/sam.cmd"
    samCmdFile = open(samCmd,"w")
    samCmdFile.write("#!/bin/csh -fe\n")
    samCmdFile.write("cd " + dir + "\n")
    samCmdFile.write(TOOLS + "/R_shell/sam.R --ac 1 -i sam.input >& sam.output\n")
    samCmdFile.write("cat sam.output | grep '.*\t.*\t.*\t.*\t' | sed -e 1d | sort >& sam.results\n")
    samCmdFile.close()
    os.chmod(samCmd, 0755);
    p = subprocess.Popen(samCmd, shell=True)
    out =  dir + "/sam.results"
    ProcessMap[p.pid] = [key, out]
    return



para = None
log = None

def lop(s):
    if s[-4] == ".":
            return s[:-4]
    return s

def main(argv):
    runNulls = False

    # initialize arguments

    if len(argv) < 3:
        die("insufficient arguments")
    
    #try:
    subprocess.check_call("which R", shell=True)
    #except Exception, e:
    #    die("Unable to find R")
    
    outfile = None
    phenotypeFile = argv[1]
    featuresFile =  argv[2]
    #BUG: need proper command line parsing
    if len(argv) > 3:
        if argv[3] == "--runNulls":
            runNulls = True
        else:
            outfile = argv[3]
    

    output = os.getcwd() + "/OCCAM" + "__" + lop(os.path.basename(phenotypeFile)) + "__" + lop(os.path.basename(featuresFile))
    return run_occam(phenotypeFile, featuresFile, output, outfile=outfile, runNulls=runNulls)
    
def run_occam(phenotypeFile, featuresFile, output, outfile=None, runNulls=False):
    random.seed(3.14159)  # first things first. ensure repeatability

    directory = output 

    global para, log;
    if not os.path.exists(directory):
        os.mkdir(directory)
    para = open( directory + "/parasol.joblist","w")
    log = open( directory + "/report.log","w")
    # log = sys.stderr

    hypotheses = directory + "/hypotheses"
    if not os.path.exists(hypotheses):
        os.mkdir(hypotheses)

    phenotype = PhenotypeTable(DichotomizedPhenotypes(phenotypeFile))
    featureTable = FeatureTable(open(featuresFile), phenotype.cellNames())
    harvester = Harvester()
    

    phenotype.runCommandForeachPhenotypeMatching(featureTable, hypotheses)
    harvester.harvestProcessMap()
    harvester.write(directory+"/results.tab")
    
    if outfile is not None:
        shutil.copyfile(directory+"/results.tab", outfile)

    if runNulls:
        nullHypotheses = directory + "/nullHypotheses"
        if not os.path.exists(nullHypotheses):
            os.mkdir(nullHypotheses)

        nc = NullCracker(featureTable)
        global ProcessMap
        global ShortCircuit
        for f in harvester.successList:
            pid = os.fork()
            if pid:
                ProcessMap[pid] = [f, nullHypotheses]
                if len(ProcessMap.keys()):
                    try:
                        (pid,status) = os.wait()
                        phe = ProcessMap[pid]
                        del ProcessMap[pid]
                    except:
                        ProcessMap = {}
            else:
                log.close()
                para.close()
                if not os.path.exists(nullHypotheses + "/" + f):
                    os.mkdir(nullHypotheses + "/" + f)
                log = open( nullHypotheses + "/" + f + "/report.log","w")
                para = open( nullHypotheses + "/" + f + "/report.log","w")
                nc.processPhenotype(f, hypotheses, nullHypotheses)
                harvester = Harvester()
                harvester.harvestProcessMap()
                harvester.write(nullHypotheses+"/"+f+".tab")
                ShortCircuit = ShortCircuit - 1
                if ShortCircuit == 0:
                   return(0)
    return(0)


class Table:
    def __init__(this, input, validate):
        this.lines = []

        # read in the lines
        lineNumber = 0
        for line in input:
           if line[-2] == "\r":
               line = line[:-2]
           else:
               line = line[:-1]
           lineAsWords = line.split("\t")
           this.lines.append(lineAsWords)
           lineNumber = lineNumber + 1
           if (lineNumber%1000) == 0:
               print >>log,"file",input.name,"line",lineNumber


        # find the row names
        this._rowMap = {}
        rowLabels = []
        for i in xrange(1,len(this.lines)):
            # rowName = this.lines[i][0].upper().replace("-","_")
            rowName = this.lines[i][0]
            if rowName == "NA":
                    rowName = "NA_"+str(i)
            rowLabels.append(rowName)
            this._rowMap[rowName] = i
        this._rowLabels = rowLabels


        # find the column names
        columnLabels = []
        this._columnMap = {}
        for column_i in  xrange(1, len(this.lines[0])):
            t = this.lines[0][column_i]
            n = t.find(" ")  # HACK to remove the column header annotations
            if n > 0:
                    t = t[:n]
            columnLabels.append(t)
            this._columnMap[t] = column_i
        this._columnLabels = columnLabels

        overlap = set(validate) & set(columnLabels) 

        if len(overlap) == 0:

            # The overlap validation test may have failed because we were handed tissue sample barcodes and patient sample barcodes.
            # remove the tissue sample barcode
            # This won't always work

            allTissue = True
            for c in columnLabels:
                if not c.endswith("-01A"):
                    allTissue = False
            if allTissue:
                print >>log,"Removing -01A tissue sample barcode to normalize features and tissue"
                for i,c in enumerate(this.lines[0]):
                    if c.endswith("-01A"):
                        this.lines[0][i] = c[:-4]
                for i,c in enumerate(columnLabels):
                    columnLabels[i] = c[:-4]

            overlap = set(validate) & set(columnLabels) 
            if len(overlap) == 0:
                die("No overlaps between phenotype and featurefile")

        this._columns = []
        for e in this.lines[0]:
            this._columns.append([])
        for line in this.lines:
           for i,e in enumerate(line):
                this._columns[i].append(e)
        print >>log,"finished loading",input.name


    def columnLabels(this):
        return this._columnLabels

    def rowLabels(this):
        return this._rowLabels

    def columns(this):
        return this._columns

    def valueAtRowColumn(this, row, column):
        try:
            return this.lines[this._rowMap[row]][this._columnMap[column]]
        except:
            return None


class FeatureTable(Table):
    def featureNames(this): # here feature names are stored in rows
       return this.rowLabels()

    def cellNames(this): # here, cell names are stored in columns
       return this.columnLabels()

    def valueAtGeneCell(this, featureName, cellName):
       value = this.valueAtRowColumn(featureName, nick(cellName)) 
       if value: return value
       return this.valueAtRowColumn(featureName, cellName)


def ZeroOne(value):
   if value == "+": return "1" 
   elif value == "-": return "0"
   else: return value

def nick(s):
   # t =  s.rsplit("-", 1)[0]
   # force exact match
   t =  s
   return t
   

class PhenotypeTable:
    def __init__(this, dich):
        this.phenotype = dich
        this.phenotypeNameToCellMap = {}
        for phenotypeName in this.phenotypeNames():
            this.phenotypeNameToCellMap[phenotypeName] = set()
            for cellName in this.cellNames():
                if this.phenotype.valueAtRowColumn(cellName, phenotypeName):
                    this.phenotypeNameToCellMap[phenotypeName].add(cellName) 


    def phenotypeNames(this): # phenotypes are stored in columns
       return this.phenotype.columnLabels()

    def cellNames(this): # here, cell names are stored in rows
       return this.phenotype.rowLabels()


    def valueAtPhenotypeTable(this, cellName, phenotypeName):
       return this.phenotype.valueAtRowColumn(cellName, phenotypeName)



    def runCommandForeachPhenotypeMatching(this, featureTable, workingDir):
      names = this.phenotypeNames()
      k = len(names)

      for i in xrange(k):
          phenotypeName = names[i]
          sys.stdout.write( "   " + str(i)+" of " + str(k)+ " " + phenotypeName+"\r")

          dir = workingDir + "/"+ phenotypeName 
          if not os.path.exists(dir):
              os.mkdir(dir)
   
          negCount = 0
          posCount = 0
          has = missing = 0
          pheSet = set(this.phenotypeNameToCellMap[phenotypeName])
          cellSet = set(featureTable.cellNames())
          hasSet = pheSet & cellSet
          thisFeaturesCellNames = list(hasSet)

          has = len(hasSet)
          misSet = pheSet - cellSet
          missing = len(misSet)

          samFilename = dir + "/sam.input"
          reportFilename = dir + "/sam.report"

          samFile = open(samFilename, "w")
          reportFile = open(reportFilename, "w")


          # TOP INFORMATION 
          # write out the header for each file
          samFile.write("\t")

          # VARIANT HEADER INFORMATION 
          tween = False
          for cellName in thisFeaturesCellNames:
              value = this.phenotype.valueAtRowColumn(cellName, phenotypeName)

              if tween:
                  samFile.write("\t")

              zeroOne = ZeroOne(value)
              if zeroOne == "1":   posCount = posCount + 1
              elif zeroOne == "0": negCount = negCount + 1

              samFile.write(zeroOne)
              reportFile.write(cellName + "\t"+ str(this.valueAtPhenotypeTable(cellName, phenotypeName)) + "\t" + zeroOne+"\n")

              tween = True

          samFile.write("\n")
          reportFile.close()

          # GENE DATA INFORMATION 
          lineNo = 1
          nValues = 0

          for featureName in featureTable.featureNames():
              samFile.write(featureName)

              for cellName in thisFeaturesCellNames:
                  samFile.write("\t")
                  value = featureTable.valueAtGeneCell(featureName, cellName)
                  if value != None:
                      nValues = nValues + 1
                      value = str(value)
                      samFile.write(value)
              samFile.write("\n")
          if nValues == 0:
             print >>log,  phenotypeName, "HAS ZERO VALUES, SKIPPING"
          else:
             print >>log,  phenotypeName,
          samFile.close()

          if True or negCount >= 5 and posCount >= 5 and dir:
              runSam(phenotypeName, dir)
              global ShortCircuit
              ShortCircuit = ShortCircuit - 1 
              if ShortCircuit == 0:
                 return
          else:
              print >>log, phenotypeName, "Insufficient samples , SKIPPING -("+str( negCount) + ") +("+str(posCount)+ ")"
               

class Harvester:
    def __init__(this):
        this.features = set()
        this.NA = set()
        this.columns = []
        this.noValue = ""
        this.successList = []
        return


    def harvestProcessMap(this):
        global ProcessMap
        while len(ProcessMap.keys()) > 0:
            try:
                (pid,status) = os.wait()
                phe = ProcessMap[pid]
                del ProcessMap[pid]
                if status == 0:
                    this.read(phe, 1)
                else:
                    print >>log, "Subprocess returned non-zero status", pid, status, phe
            except:
                # clean up any residue
                for phe in ProcessMap.values():
                    this.read(phe, 1)
                ProcessMap = {}

    def read(this, phe, columnNumber):
        " read in filename and extract the zero'th colum and columnNumber into our data"
        try:
            print >>log, phe[0], 
            lineNumber = 0
            (header, filename) = phe
            column = {}
            for line in open(filename):
               lineNumber = lineNumber + 1
               line = line[:-1]
               words = line.split("\t")
               if len(words) >  1:
                   feature = words[0]
                   this.features.add(feature)
                   value =  words[columnNumber]
                   if value == "NA":
                       this.NA.add(feature)
                   elif float(value) > 20:
                       pdb.set_trace()
                   else:
                       column[feature] = value
            if lineNumber < 10 or len(column.keys()) < 10:
                raise Exception("insufficient harvest")

            this.columns.append(column)
            this.successList.append(header)
        except:
            print >>sys.stdout,"FAILURE: SAM processing for ",header

            print >>log,"FAILURE: SAM processing for ",phe
            print >>log,"BEGIN"
            print >>log, open(filename.replace("sam.results","sam.output")).read()
            print >>log,"END"
        """
        except as inst:
            print >>sys.stderr,"FAILURE: SAM processing for ",header

            print >>log,"FAILURE: SAM processing for ",phe
            print >>log, type(inst)     # the exception instance
            print >>log, inst           # __str__ allows args to printed directly
            print >>log,"BEGIN"
            print >>log, open(filename.replace("sam.results","sam.output")).read()
            print >>log,"END"
        """
    
    def write(this, filename):
         " write out the table as a tab delimited file"
         print >>log,"\n results are in\n ", filename
         file = open(filename,"w")

         file.write("id")
         for header in this.successList:
            file.write("\t")
            file.write(header)
         file.write("\n")

         for feature in sorted(this.features):
            file.write(feature)
            for c in this.columns:
                file.write("\t")
                if feature in c:
                    file.write(c[feature])
            file.write("\n")
         file.close()



if __name__ == "__main__":
        #import cProfile
        #cProfile.run("main(sys.argv)")
        main(sys.argv)

