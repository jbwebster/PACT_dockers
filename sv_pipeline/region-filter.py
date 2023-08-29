#!/usr/bin/env python

#########
# Apply region based filters to SURVIVOR output
#
########

import sys
import re
import argparse

class Region:
	def __init__(self, chrm, start, end):
		if "chr" not in chrm:
			chrm = "chr" + str(chrm)
		self.chrm = chrm
		self.start = start
		self.end = end

class Variant:
	def __init__(self, chrma, starta, chrmb, startb, line):
		if "chr" not in chrma:
			chrma = "chr" + str(chrma)
		self.chrma = chrma
		self.starta = int(starta)
		if "chr" not in chrmb:
			chrmb = "chr" + str(chrmb)
		self.chrmb = chrmb
		self.startb = int(startb)
		self.line = line
		self.keep = True

	def atLeastOne(self, chrmalist, chrmblist):
		if not chrmalist and not chrmblist:
			return False
		if chrmalist:
			for region in chrmalist:
				if self.starta >= region.start and self.starta <= region.end:
					return True
		if chrmblist:
			for region in chrmblist:
				if self.startb >= region.start and self.startb <= region.end:
					return True
		return False

	def noMoreThanOne(self, chrmalist, chrmblist):
		if not chrmalist or not chrmblist:
			return True
		matches = 0
		if chrmalist:
			for region in chrmalist:
				if self.starta >= region.start and self.starta <= region.end:
					matches += 1
		if chrmblist:
			for region in chrmblist:
				if self.startb >= region.start and self.startb <= region.end:
					matches += 1
		if matches > 1:
			return False
		return True

	def noOverlap(self, chrmalist, chrmblist):
		if not chrmalist and not chrmblist:
			return True
		if chrmalist:
			for region in chrmalist:
				if self.starta >= region.start and self.starta <= region.end:
					return False
		if chrmblist:
			for region in chrmblist:
				if self.startb >= region.start and self.startb <= region.end:
					return False
		return True

def createParser():
	parser = argparse.ArgumentParser(description="Region filter")
	parser.add_argument("-v", "--vcf", action="store", dest="vcf", help="SURVIVOR vcf", required=True)
	parser.add_argument("-t", "--target", action="store", dest="target", help="Target bed", required=True)
	parser.add_argument("-l", "--lowc", action="store", dest="lowc", help="Low complexity bed")
	parser.add_argument("-b", "--blacklist", action="store", dest="blacklist", help="Blacklist bed")
	return parser

def readRegions(bed, isTarget):
	regions = {}
	isHeader = True
	with open(bed, "r") as f:
		for line in f:
			if isHeader:
				try:
					x = int(line.split("\t")[1])
					isHeader = False
				except:
					isHeader = True
				
			if not isHeader:
				line = line.strip()
				fields = line.split("\t")
				mod = 0
				if isTarget:
					mod = 200
				region = Region(fields[0], int(fields[1]) - mod, int(fields[2]) + mod)
				if region.chrm in regions.keys():
					l = regions[region.chrm]
					l.append(region)
					regions[region.chrm] = l
				else:
					l = [region]
					regions[region.chrm] = l
	f.close()	
	return regions

def parseInfo(info,line):
	chrma = info.split("_")
	if len(chrma) > 2:
		chrma = chrma[0] + chrma[1]
	else:
		chrma = chrma[0]
	starta = int(info.split("-")[0].split("_")[-1])
	chrmb = info.split("-")[1].split("_")
	if len(chrmb) > 2:
		chrmb = chrmb[0] + chrmb[1]
	else:
		chrmb = chrmb[0]
	startb = int(info.split("-")[1].split("_")[-1])
	variant = Variant(chrma,starta,chrmb,startb,line)	
	return variant

def readVCF(vcf):
	headers = []
	variants = []
	with open(vcf, "r") as f:
		for line in f:
			line = line.strip()
			if line[0] == "#":
				headers.append(line.strip())
			else:
				fields = line.split("\t")
				toola = fields[9]
				toolb = fields[10]
				toolc = fields[11]
				toolainfo = toola.split(":")[-1]
				toolbinfo = toolb.split(":")[-1]
				toolcinfo = toolc.split(":")[-1]
				chrma = None
				starta = None
				chrmb = None
				startb = None
				added = False
				if "NA" not in toolainfo:
					if "," in toolainfo:
						toolainfo = toolainfo.split(",")[-1]
					variant = parseInfo(toolainfo, line)
					variants.append(variant)
					added = True
				if "NA" not in toolbinfo and not added:
					if "," in toolbinfo:
						toolbinfo = toolbinfo.split(",")[-1]
					variant = parseInfo(toolbinfo, line)
					variants.append(variant)
					added = True
				if "NA" not in toolcinfo and not added:
					if "," in toolcinfo:
						toolcinfo = toolcinfo.split(",")[-1]
					variant = parseInfo(toolcinfo, line)
					variants.append(variant)
	f.close()
	return headers, variants

def atLeastOne(variants, regions):
	for variant in variants:
		ra = False
		rb = False
		if variant.chrma in regions.keys():
			ra = regions[variant.chrma]
		if variant.chrmb in regions.keys():
			rb = regions[variant.chrmb]	
		if not variant.atLeastOne(ra,rb):
			variant.keep = False
	return variants		

def noMoreThanOne(variants, regions):
	for variant in variants:
		ra = False
		rb = False
		if variant.chrma in regions.keys() and variant.keep:
			ra = regions[variant.chrma]
		if variant.chrmb in regions.keys() and variant.keep:
			rb = regions[variant.chrmb]
		if not variant.noMoreThanOne(ra,rb):
			variant.keep = False
	return variants

def noOverlap(variants, regions):
	for variant in variants:
		ra = False
		rb = False
		if variant.chrma in regions.keys() and regions.keep:
			ra = regions[variant.chrma]
		if variant.chrmb in regions.keys() and variant.keep:
			rb = regions[variant.chrmb]
		if not variant.noOverlap(ra,rb):
			variant.keep = False
	return variants
	
def printOutput(headers,variants):
	for header in headers:
		print(header.strip())
	for variant in variants:
		if variant.keep:
			print(variant.line.strip())


def main(args):
	parser = createParser()	
	args = parser.parse_args()
	targets = readRegions(args.target, True)
	headers, variants = readVCF(args.vcf)
	variants = atLeastOne(variants, targets)
	if args.lowc:
		lowc = readRegions(args.lowc, False)
		variants = noMoreThanOne(variants, lowc)
	if args.blacklist:
		blacklist = readRegions(args.blacklist, False)
		variants = noOverlap(variants, blacklist)
	printOutput(headers,variants)
	

if __name__ == '__main__':
	sys.exit(main(sys.argv))
