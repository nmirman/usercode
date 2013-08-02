#!/usr/bin/env python

import sys
import optparse
import re
from FWCore.PythonUtilities.LumiList import LumiList

if __name__ == '__main__':

   parser = optparse.OptionParser ("Usage: %prog alpha.json")
   parser.add_option ('--sec', dest='sec', type='int', default=0, help='number of sections')
   parser.add_option ('--outbase', dest='outbase', type='string', default = 'temp',
         help='base for output filenames')
   parser.add_option ('--outend', dest='outend', type='string', default = [],
         help='ending for output filenames')

   # required parameters
   (options, args) = parser.parse_args()
   if len (args) != 1:
      raise RuntimeError, "Must provide exactly one input file"

   if options.sec < 2:
      raise RuntimeError, "Need at least 2 sections."

   commaRE = re.compile (r',')
   
   alphaList = LumiList (filename = args[0]) # Read in first JSON file
   allLumis = alphaList.getLumis()

   count_lumis = 0
   for (run, lumi) in allLumis:
      count_lumis += 1

   lumis_per_sec = count_lumis/options.sec
   print 'Found %s lumis in file.' % count_lumis
   print 'Splitting into %s sections with %s lumis each.' % (options.sec, lumis_per_sec)

   count_lumis = 0
   for sec in range(options.sec):
      
      ibegin = sec*lumis_per_sec
      iend = ibegin + lumis_per_sec
      if sec == options.sec - 1:
         iend = len(allLumis)

      tempList = LumiList (lumis = allLumis[ibegin:iend])
      print 'Part %s: num lumis = %s' % (sec+1, len(tempList.getLumis()))
      count_lumis += len(tempList.getLumis())

      filename = '%s_%s_%s' % (options.outbase, sec+1, options.outend)
      tempList.writeJSON (filename)

print 'Total lumis = %s' % count_lumis
