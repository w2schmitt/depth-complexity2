#!/usr/bin/python

import sys, os, getopt, io, re, json, operator

def main(argv):
	inputDir = "."

	# parse input
	try:
		opts, args = getopt.getopt(argv, "hi:",["idir="])
	except getopt.GetoptError:
		print 'AvgResultsInfo.py -i <inputDir>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'AvgResultsInfo.py -i <inputDir>'
			sys.exit()
		elif opt in ("-i", "--idir"):
			inputDir = arg

	# get list of all files
	listFiles = []
	for root, dir, files in os.walk(inputDir):
		for f in files:
			if (f == "info.txt"):
				listFiles.append( os.path.join(root, f)	)
	listFiles.sort()
	
	# parse
	n = float(len(listFiles))
	acumulator = {"s" : 0, "ms" : 0, "tri" : 0, "vert" : 0}
	maxEntry = {"s" : 0, "ms" : 0, "tri" : 0, "vert" : 0}
	minEntry = {"s" : 0, "ms" : 99999, "tri" : 0, "vert" : 0}	
	for path in listFiles: 
		currentEntry = {}
		f = open(path, 'r')
		text = filter(None, re.split(' |\n', f.read())[4:])
		for i in xrange(0,len(text)-1,2):
			field = text[i]
			value = float(text[i+1])
			acumulator[field] += float(value)
			currentEntry[text[i]] = float(text[i+1])
		if "s" in currentEntry:
			currentEntry["ms"] += currentEntry["s"]*1000
			del currentEntry["s"]
		if (currentEntry["ms"] > maxEntry["ms"]):
			maxEntry = currentEntry
		if (currentEntry["ms"] < minEntry["ms"]):
			minEntry = currentEntry
		
	acumulator["ms"] += acumulator["s"]*1000
	del acumulator["s"]

	print minEntry
	print maxEntry
	print acumulator

	# max value	
	acumulator.update((x, float(y)/n) for x, y in acumulator.items())

	print "\nAverage Values:"
	print json.dumps(acumulator, indent=4)

# ----------------------------- main
if __name__ == "__main__":
	main(sys.argv[1:])