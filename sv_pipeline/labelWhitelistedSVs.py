#!/usr/bin/python

import sys
import argparse

def createParser():
	parser = argparse.ArgumentParser(description="Label whitelisted SVs")
	parser.add_argument("-w", action="store", dest="whitelist")
	parser.add_argument("-b", action="store", dest="bedpe")
	return parser

def readWhitelist(infile):
	whitelist = []	
	with open(infile, "r") as f:
		for line in f:
			fields = line.split("\t")
			if len(fields) == 4:
				classes = ["DUP","INV","DEL","BND"]
				if fields[3] in classes:
					hotspot = [fields[0].strip(), int(fields[1].strip()), int(fields[2].strip()), fields[3].strip()]
				else:
					hotspot = [fields[0].strip(), int(fields[1].strip()), int(fields[2].strip()), fields[3].strip()]
				whitelist.append(hotspot)
			if len(fields) == 3:
				hotspot = [fields[0].strip(), int(fields[1].strip()), int(fields[2].strip())]
				whitelist.append(hotspot)
	f.close()
	return whitelist	

def applyWhitelist(whitelist, bedpe):
	with open(bedpe, "r") as f:
		for line in f:
			if "/" in line:
				line = line.strip()
				fields = line.split("\t")
				pos1 = (int(fields[1]) + int(fields[2])) / 2
				pos2 = (int(fields[4]) + int(fields[5])) / 2
				chr1 = fields[0]
				chr2 = fields[3]
				svtype = fields[6].split("/")[1]
				isWhitelisted = False
				for hotspot in whitelist:
					if chr1 == hotspot[0] and pos1 >= hotspot[1] and pos1 <= hotspot[2]:
						if len(hotspot) == 3 or (len(hotspot) == 4 and svtype == hotspot[3]):
							isWhitelisted = True
					if chr2 == hotspot[0] and pos2 >= hotspot[1] and pos2 <= hotspot[2]:
						if len(hotspot) == 3 or (len(hotspot) == 4 and svtype == hotspot[3]):
							isWhitelisted = True
				info = fields[20]
				if isWhitelisted:
					info = "WHITELISTED=TRUE;" + info
				else:
					info = "WHITELISTED=FALSE;" + info
				fields[20] = info
				out = "\t".join(fields)
				print(out)	
	f.close()

def main(args):
	parser = createParser()
	args = parser.parse_args()
	
	bedpe = args.bedpe	
	if args.whitelist is not None:
		whitelist = readWhitelist(args.whitelist)
		applyWhitelist(whitelist, bedpe)
	else: # If no whitelist is supplied, just output the original bedpe file
		with open(bedpe, "r") as f:
			for line in f:
				print(line.strip())
		f.close()

if __name__ == '__main__':
	sys.exit(main(sys.argv))



