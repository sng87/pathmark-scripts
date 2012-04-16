#!/usr/bin/env python
"""tab2html.py: 

Usage:
  tab2html.py [options] tabFile htmlFile

Options:
  -q            run quietly
"""
import os, os.path, sys, getopt, re
import mData

verbose = True

htmlHead = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
<html>

    <head>
        <title>HTML Table</title>
        <script type="text/javascript" src="../js/jquery-1.6.1.min.js"></script>
        <script type="text/javascript" src="../js/jquery.tablesorter.js"></script>
        <script type="text/javascript">
            $(document).ready(function() 
                { 
                    $("#htmlTable").tablesorter(); 
                }
            );
        </script>
    </head>
    <body>
"""

htmlIndexItem = """        <a href="%s">%s</a><br>
"""

htmlCategory = """    <h1>%s</h1>
"""

htmlTableHead = """    <table id="htmlTable" class="tablesorter" border="1">
        <thead>
        <tr>
            <th>%s</th>
        </tr>
        </thead>
        <tbody>
"""

htmlTableItem = """        <tr>
            <td>%s</td>
            <td>%s</td>
        </tr>
"""

htmlTableTail = """        </tbody>
    </table>
"""

htmlTail = """    </body>
</html>
"""

def usage(code = 0):
    print __doc__
    if code != None:
        sys.exit(code)

def log(msg):
    if (verbose):
        sys.stderr.write(msg)

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "i:q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if (len(args) != 1) & (len(args) != 2):
        print "incorrect number of arguments"
        usage(1)
    
    if (len(args) == 1):
        tabFile = sys.stdin
        htmlFile = args[0]
    else:
        tabFile = args[0]
        htmlFile = args[1]
    
    global verbose
    for o, a in opts:
        if o == "-q":
            verbose = False
    
    ## read tabFile
    (tabData, tabCols, tabRows) = mData.rCRSData(tabFile, retFeatures = True)
    
    ## write htmlFile
    f = open("%s" % (htmlFile), "w")
    f.write(htmlHead)
    tabCols = ["id"] + tabCols
    f.write(htmlTableHead % ("</th>\n            <th>".join(tabCols)))    
    for i in tabRows:
        tabList = []
        for j in tabCols[1:]:
            tabList.append(str(tabData[j][i]))
        f.write(htmlTableItem % (i, "</td>\n            <td>".join(tabList)))
    f.write(htmlTableTail)
    f.write(htmlTail)
    f.close()

if __name__ == "__main__":
    main(sys.argv[1:])