import math, os, os.path, sys, getopt, re

class Pathway:
    def __init__(self, nodes, interactions):
        self.nodes = nodes
        self.interactions = interactions
    def removeNode(self, node):
        (self.nodes, self.interactions) = removeNode(node, self.nodes, self.interactions)

def log(msg, die = False):
    sys.stderr.write(msg)
    if (die):
        sys.exit(1)

def openAnyFile(inf):
    """performs an open() on a file or url"""
    if inf.startswith("http"):
        import urllib2
        stream = urllib2.urlopen(inf)
    else:
        stream = open(inf, 'r')
    return stream

def retColumns(inf, delim = "\t"):
    """returns the columns of a .tsv"""
    f = openAnyFile(inf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank header\n", die = True)
    line = line.rstrip("\r\n")
    return(re.split(delim, line)[1:])

def rCRSData(inf, appendData = dict(), delim = "\t", retFeatures = False, debug = False):
    """reads .tsv into a [col][row] dictionary"""
    inData = dict()
    colFeatures = []
    rowFeatures = []
    ## copy appendData
    for col in appendData.keys():
        inData[col] = dict()
        for row in appendData[col].keys():
            inData[col][row] = [appendData[col][row]]
    ## read header
    f = openAnyFile(inf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    if debug:
        log("%s\nLENGTH: %s\n" % (line, len(pline)))
    colFeatures = pline[1:]
    for col in colFeatures:
        if col not in inData:
            inData[col] = dict()
    ## read data
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        rowFeatures.append(pline[0])
        if debug:
            log("%s\nLENGTH: %s\n" % (line, len(pline)))
        if len(pline) != (1+len(colFeatures)):
            log("ERROR: length of line does not match the rest of the file\n", die = True)
        for i, col in enumerate(colFeatures):
            row = pline[0]
            if row not in inData[col]:
                inData[col][row] = []
            if pline[i+1] == "":
                inData[col][row].append("NA")
            else:            
                inData[col][row].append(pline[i+1])
    f.close()
    ## average entries
    for col in inData.keys():
        for row in inData[col].keys():
            inData[col][row] = mean(inData[col][row], null = inData[col][row][0])
    if retFeatures:
        return(inData, colFeatures, rowFeatures)
    else:
        return(inData)

def rwCRSData(outf, inf, delim = "\t", null = "NA", useCols = None, useRows = None, colMap = {}, rowMap = {}, rcMap = None):
    """reads and writes .tsv"""
    colFeatures = []
    ## read header
    f = openAnyFile(inf)
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    lineLength = len(pline)
    colIndex = {}
    for i, col in enumerate(pline[1:]):
        colIndex[col] = i
        if useCols != None:
            if col not in useCols:
                continue
        colFeatures.append(col)
    ## write header
    if os.path.exists(outf):
        o = open(outf, "a")
    else:
        o = open(outf, "w")
        o.write("id%s\n" % (delim+delim.join(colFeatures)))
    ## read and write data
    rowCount = 0
    for line in f:
        if line.isspace():
            continue
        rowCount += 1
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        if pline[0] in rowMap:
            mrow = rowMap[pline[0]]
        else:
            mrow = pline[0]
        if useRows != None:
            if mrow not in useRows:
                continue
        if len(pline) != lineLength:
            log("ERROR: length of line does not match the rest of the file\n", die = True)
        else:
            o.write("%s" % (mrow))
        if rcMap is not None:
            colMap = rcMap[mrow]
        for col in colFeatures:
            if col in colMap:
                mcol = colMap[col]
            else:
                mcol = col
            if pline[colIndex[mcol]+1] == "":
                o.write("%s" % (delim+null))
            else:            
                o.write("%s" % (delim+pline[colIndex[mcol]+1]))
        o.write("\n")
    f.close()
    o.close()

def rList(inf, header = False):
    """read 1 column list"""
    inList = []
    f = openAnyFile(inf)
    if header:
        f.readline()
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\t\r\n")
        inList.append(line)
    f.close()
    return(inList)

def getListIndices(inItem, inList):
    """returns indices of the occurence of inItem in inList"""
    indices = []
    for i, item in enumerate(inList):
        if item == inItem:
            indices.append(i)
    return(indices)

def floatList(inList):
    """returns only numeric elements of a list"""
    outList = []
    for i in inList:
        try:
            fval = float(i)
            if fval != fval:
                raise ValueError
            outList.append(fval)
        except ValueError:
            continue
    return(outList)

def mean(inList, null = "NA"):
    """Calculates mean"""
    fList = floatList(inList)
    if len(fList) == 0:
        mean = null
    else:
        mean = sum(fList)/len(fList)
    return (mean)

def mean_std(inList, sample = True):
    """Calculates mean and std"""
    
    cList = floatList(inList)
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

def addLink(a, b, pNodes, pInteractions, gNodes, gInteractions):
    if a not in pNodes:
        pNodes[a] = gNodes[a]
    if b not in pNodes:
        pNodes[b] = gNodes[b]
    if a not in pInteractions:
        pInteractions[a] = dict()
    pInteractions[a][b] = gInteractions[a][b]
    return(pNodes, pInteractions)

def rPathway(inf, reverse = False, retProteins = False, delim = "\t"):
    """read UCSC pathway tab"""
    proteins = set()
    readPathway = Pathway(dict(), dict())
    f = open(inf, "r")
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        if len(pline) == 2:
            readPathway.nodes[pline[1]] = pline[0]
            if pline[0] == "protein":
                proteins.update([pline[1]])
        elif len(pline) == 3:
            if reverse:
                if pline[1] not in readPathway.interactions:
                    readPathway.interactions[pline[1]] = dict()
                if pline[0] not in readPathway.interactions[pline[1]]:
                    readPathway.interactions[pline[1]][pline[0]] = pline[2]
                else:
                    readPathway.interactions[pline[1]][pline[0]] += ";"+pline[2]
            else:
                if pline[0] not in readPathway.interactions:
                    readPathway.interactions[pline[0]] = dict()
                if pline[1] not in readPathway.interactions[pline[0]]:
                    readPathway.interactions[pline[0]][pline[1]] = pline[2]
                else:
                    readPathway.interactions[pline[0]][pline[1]] += ";"+pline[2]
        else:
            print >> sys.stderr, "ERROR: line length not 2 or 3: \"%s\"" % (line)
            sys.exit(1)
    f.close()
    if retProteins:
        return(readPathway.nodes, readPathway.interactions, proteins)
    else:
        return(readPathway.nodes, readPathway.interactions)

def wPathway(outf, outNodes, outInteractions, useNodes = None):
    """write UCSC pathway.tab"""
    f = open(outf, "w")
    if useNodes == None:
        useNodes = outNodes.keys()
    for i in useNodes:
        if i not in outNodes:
            continue
        f.write("%s\t%s\n" % (outNodes[i], i))
    for i in useNodes:
        if i not in outInteractions:
            continue
        for j in outInteractions[i].keys():
            if j not in useNodes:
                continue
            for k in re.split(";", outInteractions[i][j]):
                f.write("%s\t%s\t%s\n" % (i, j, k))
    f.close()

def wSIF(writeFile, writeInteractions, useNodes = None):
    """write .sif"""
    f = open(writeFile, "w")
    if useNodes == None:
        for i in writeInteractions.keys():
            for j in writeInteractions[i].keys():
                for k in re.split(";", writeInteractions[i][j]):
                    f.write("%s\t%s\t%s\n" % (i, k, j))
    else:
        for i in useNodes:
            if i not in writeInteractions:
                continue
            for j in writeInteractions[i].keys():
                if j not in useNodes:
                    continue
                for k in re.split(";", writeInteractions[i][j]):
                    f.write("%s\t%s\t%s\n" % (i, k, j))
    f.close()

def getComponentMap(pNodes, pInteractions):
    """create the dictionary componentMap from interaction map"""
    rpInteractions = reverseInteractions(pInteractions)
    componentMap = dict()
    for i in pNodes.keys():
        if pNodes[i] != "complex":
            continue
        componentMap[i] = []
        if i not in rpInteractions:
            continue
        for j in rpInteractions[i]:
            if rpInteractions[i][j] == "component>":
                componentMap[i].append(j)
    return(componentMap)

def filterComplexesByGeneSupport(allNodes, forInteractions, revInteractions, typeMap, componentMap, threshold = 0.5):
    """remove complexes by percent support"""
    ## mark complexes in allNodes as unvisitedComplexes
    unvisitedComplexes = set()
    for i in allNodes:
        if allNodes[i] == "complex":
            unvisitedComplexes.update([i])
    def keepMajority(complex, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions, typeMap, componentMap, threshold = 0.5):
        unvisitedComplexes.discard(complex)
        complexSubunits = []
        otherSubunits = []
        for i in componentMap[complex]:
            if typeMap[i] == "complex":
                complexSubunits.append(i)
            else:
                otherSubunits.append(i)
        n = len(complexSubunits)+len(otherSubunits)
        complexSubunitsInNet = []
        otherSubunitsInNet = []
        for i in complexSubunits:
            if i in allNodes:
                if i in recursedComplexes:
                    indices = getListIndices(i, recursedComplexes)
                    log("WARNING: complex loop %s\n" % (recursedComplexes[indices[-1]:]))
                    continue
                recursedComplexes.append(i)
                (logical, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions) = keepMajority(i, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions, typeMap, componentMap, threshold = threshold)
                if logical:
                    complexSubunitsInNet.append(i)
        for i in otherSubunits:
            if i in allNodes:
                otherSubunitsInNet.append(i)
        m = len(complexSubunitsInNet)+len(otherSubunitsInNet)
        if m <= threshold*n:
            (allNodes, forInteractions) = removeNode(complex, allNodes, forInteractions)
            return(False, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions)
        else:
            return(True, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions)

    ## visit complexes while there are still unvisitedComplexes
    while len(unvisitedComplexes) > 0:
        complex = list(unvisitedComplexes)[0]
        recursedComplexes = [complex]
        (logical, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions) = keepMajority(complex, unvisitedComplexes, recursedComplexes, allNodes, forInteractions, revInteractions, typeMap, componentMap, threshold = threshold)
    return(allNodes, forInteractions)

def reverseInteractions(pInteractions):
    """reverse interaction mapping"""
    rpInteractions = dict()
    for i in pInteractions.keys():
        for j in pInteractions[i].keys():
            if j not in rpInteractions:
                rpInteractions[j] = dict()
            rpInteractions[j][i] = pInteractions[i][j]
    return(rpInteractions)

def constructInteractions(nodeList, refNodes, refInteractions):
    """select concepts from list and construct Pathway"""
    outPathway = Pathway({}, {})
    for i in nodeList:
        outPathway.nodes[i] = refNodes[i]
        if i in refInteractions:
            for j in refInteractions[i].keys():
                if j in nodeList:
                    if i not in outPathway.interactions:
                        outPathway.interactions[i] = dict()
                    outPathway.interactions[i][j] = refInteractions[i][j]
    return(outPathway.nodes, outPathway.interactions)

def sortConnected(allNodes, forInteractions, revInteractions, method = "size", addData = None):
    index = 1
    mapNets = dict()
    sortedNets = []
    seenNodes = set()
    for i in allNodes.keys():
        if i in seenNodes:
            continue
        borderNodes = [i]
        currentNet = [i]
        while len(borderNodes) > 0:
            if borderNodes[0] in revInteractions:
                for j in revInteractions[borderNodes[0]].keys():
                    if j not in seenNodes:
                        seenNodes.update([j])
                        borderNodes.append(j)
                        currentNet.append(j)
            if borderNodes[0] in forInteractions:
                for j in forInteractions[borderNodes[0]].keys():
                    if j not in seenNodes:
                        seenNodes.update([j])
                        borderNodes.append(j)
                        currentNet.append(j)
            borderNodes.pop(0)
        if ("__DISCONNECTED__" not in currentNet):
            mapNets[index] = deepcopy(currentNet)
            index += 1
    indexList = mapNets.keys()
    netScore = dict()
    for i in indexList:
        if method == "size":
            netScore[i] = len(mapNets[i])
        elif method == "average":
            values = []
            for j in mapNets[i]:
                if j in addData:
                    if addData[j] != "NA":
                        values.append(abs(addData[j]))
            if len(values) > 0:
                netScore[i] = sum(values)/len(values)
            else:
                netScore[i] = 0.0
        elif method == "overlap":
            netScore[i] = len(set(mapNets[i]) & addData)
    indexList.sort(lambda x, y: cmp(netScore[y], netScore[x]))
    for i in indexList:
        sortedNets.append(mapNets[i])
    return(sortedNets)

def removeNode(node, pNodes, pInteractions):
    """remove a node and its interactions from a Pathway"""
    rpInteractions = reverseInteractions(pInteractions)
    del pNodes[node]
    if node in pInteractions:
        for element in pInteractions[node].keys():
            del pInteractions[node][element]
            if len(pInteractions[node].keys()) == 0:
                del pInteractions[node]
            del rpInteractions[element][node]
            if len(rpInteractions[element].keys()) == 0:
                del rpInteractions[element]
    if node in rpInteractions:
        for element in rpInteractions[node].keys():
            del pInteractions[element][node]
            if len(pInteractions[element].keys()) == 0:
                del pInteractions[element]
            del rpInteractions[node][element]
            if len(rpInteractions[node].keys()) == 0:
                del rpInteractions[node]
    return(pNodes, pInteractions)

def wNodeAttributes(feature, gNodes, scoreData):
    if not os.path.exists("TYPE.NA"):
        typef = open("TYPE.NA", "w")
        typef.write("TYPE (class=java.lang.String)\n")
        for node in gNodes.keys():
            typef.write("%s = %s\n" % (node, gNodes[node]))
        typef.close()
    if not os.path.exists("LABEL.NA"):
        labelf = open("LABEL.NA", "w")
        labelf.write("LABEL (class=java.lang.String)\n")
        for node in gNodes.keys():
            if gNodes[node] == "protein":
                labelf.write("%s = %s\n" % (node, node))
            else:
                labelf.write("%s = %s\n" % (node, ""))
        labelf.close()
    scoref = open(feature+"_SCORE.NA", "w")
    scoref.write("SCORE (class=java.lang.Float)\n")
    for node in gNodes.keys():
        if node in scoreData[feature]:
            if scoreData[feature][node] == "NA":
                scoref.write("%s = %s\n" % (node, "0"))
            else:
                scoref.write("%s = %s\n" % (node, scoreData[feature][node]))
        else:
            scoref.write("%s = %s\n" % (node, "0"))
    scoref.close()

def selectLink(feature, source, target, sData, pStats, filterBounds, selectionRule = "OR"):
    linkScore = []
    srcScore = []
    trgScore = []
    for i in range(len(sData.keys())):
        linkScore.append([sData[i][feature][source], sData[i][feature][target]])
    for i in range(len(sData.keys())):
        if linkScore[i][0] > pStats[i][0]+filterBounds[1]*pStats[i][1]:
            srcScore.append(2)
        elif linkScore[i][0] > pStats[i][0]+filterBounds[0]*pStats[i][1]:
            srcScore.append(1)
        else:
            srcScore.append(0)
        if linkScore[i][1] > pStats[i][0]+filterBounds[1]*pStats[i][1]:
            trgScore.append(2)
        elif linkScore[i][1] > pStats[i][0]+filterBounds[0]*pStats[i][1]:
            trgScore.append(1)
        else:
            trgScore.append(0)
    
    ## selection rule
    if selectionRule == "OR":
        if max(srcScore)+max(trgScore) >= 3:
            return(True)
    elif selectionRule == "AND":
        votes = 0
        for i in range(len(sData.keys())):
            if srcScore[i]+trgScore[i] >= 3:
                votes += 0
        if votes == len(sData.keys()):
            return(True)
    elif selectionRule == "MAIN":
        if srcScore[0]+trgScore[0] >= 3:
            return(True)
    return(False)
