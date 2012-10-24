#!/usr/bin/env python

##########################################################################
# Filename: assemblyGraphToDot.py
# Author: Lee Mendelowitz
# Date: 8/28
#
# Description:
# Parse an SGA ASQG file and output a dot file

# Usage:
# ./assemblyGraphToDot_ver2.py [-o OUTPUT_PATH] ASQG_PATH
##########################################################################
import os
import sys
import re
import argparse

parser = argparse.ArgumentParser(description='Covert SGA ASQG to dot format')
parser.add_argument('ASQGFile', help = 'ASQG File Path')
parser.add_argument('-o', '--output', help='Output File Path. Default: STDOUT')

edgeFields = [('contig1', str),
              ('contig2', str),
              ('s1', int), 
              ('e1', int),
              ('l1', int),
              ('s2', int), 
              ('e2', int),
              ('l2', int),
              ('rc', bool),
              ('numDiff', int)]
edgeFieldNames, edgeFieldTypes = zip(*edgeFields)
numEdgeFields = len(edgeFields)

##########################################
# Make an edge from the fields
# The edge is stored as a dictionary object
# Make overlap indices python friendly (by adding 1 to end index of overlap)
def makeEdge(fields):
    if len(fields) != numEdgeFields:
        raise RuntimeError('Incorrect number of fields in edge record!')
    edgeDict = dict((k,t(v)) for k,t,v in zip(edgeFieldNames, edgeFieldTypes, fields))  
    edgeDict['e1'] += 1
    edgeDict['e2'] += 1
    edgeDict['end1'] = 'B' if edgeDict['s1']==0 else 'E'
    edgeDict['end2'] = 'B' if edgeDict['s2']==0 else 'E'
    return edgeDict

##########################################
# Read an asqgFile
# Return:
#   contigDict (keys are contigIds, values is dictionary of contig attributes)
#   edgeList (list of edge dictionaries)
def parseASQG(asqgFile):
    fin = open(asqgFile)
    contigDict = {}
    edgeList = []
    for line in fin:
        fields = line.split()
        recTag = fields[0]
        if recTag=='VT':
            if len(fields) != 3:
                raise RuntimeError('Expected 3 fields in contig record: ' + str(fields))
            contigId = fields[1]
            seq = fields[2]
            assert (contigId not in contigDict)
            contigDict[contigId] = dict( [('len'  , len(seq)),
                                          ('label', contigId) ] )
        elif recTag=='ED':
            edgeList.append(makeEdge(fields[1:]))
    fin.close()
    sys.stderr.write('Read %i nodes, %i edges\n'%(len(contigDict), len(edgeList)) )
    return contigDict, edgeList

##########################################
# Write the graph specified by nodeDict and edgeList to dot file
def writeDot(nodeDict, edgeList, fout = sys.stdout):
    fout.write('digraph {\n')

    # Write the nodes with their relevent attributes: fillcolor, label
    convNode = lambda s: re.sub(r'\W','',s) # remove non-alphanumeric chars
    for node, nattr in nodeDict.iteritems():
        attrString = ''
        attrString += ' style=filled fillcolor=%s '%nattr['fillcolor'] if 'fillcolor' in nattr else ''
        attrString += ' label="%s" ' % nattr['label']
        fout.write('%s [%s];\n'%(convNode(node), attrString))

    # Write the edges wtih their relevant attributes: 
    for eattr in edgeList:
        n1, n2 = eattr['contig1'], eattr['contig2']
        attrString = ' dir=both '
        # If end is B, arrow points in (normal). If end is E, arrow points out (inv)
        if eattr['end1']=='B':
            arrowtail = 'normal'; tailport = 'w'
        else:
            arrowtail = 'inv'; tailport = 'e'
        attrString += ' arrowtail=%s tailport=%s'%(arrowtail, tailport)
        if eattr['end2']=='B':
            arrowhead = 'normal'; headport = 'w'
        else:
            arrowhead = 'inv'; headport = 'e'
        attrString += ' arrowhead=%s headport=%s '%(arrowhead, headport)
        fout.write('%s->%s [%s];\n'%(convNode(n1), convNode(n2), attrString))
    fout.write('}')

##########################################
def graphToDot(asqgFileName, fout = sys.stdout):
    contigDict, edgeList = parseASQG(asqgFileName)
    writeDot(contigDict, edgeList, fout)

##########################################
def main():
    args = parser.parse_args()

    asqgFileName = args.ASQGFile
    outputFileName = args.output

    sys.stderr.write('Inputs:\n')
    sys.stderr.write('ASQG File: %s\n'%asqgFileName)
    sys.stderr.write('Output File: %s\n'% (outputFileName if outputFileName else 'STDOUT'))

    if not os.path.exists(asqgFileName):
        raise RuntimeError('Could not find ASQG file %s'%asqgFileName)

    if outputFileName:
        fout = open(outputFileName, 'w')
    else:
        fout = sys.stdout
    graphToDot(asqgFileName, fout)

##########################################
if __name__ == '__main__':
    main()
