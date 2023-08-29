#!/usr/bin/python

import sys
import re
import argparse
import math

"""
Enforce min vaf and min coverage filters, since some make it pass other filtering steps for some reason.
Flip all SNVs from the whitelist/curated list to 'PASS' if desired
"""

def printOriginal(vcf):
	with open(vcf, "r") as f:
		for line in f:
			print(line.rstrip())	
	f.close()

def printHeader(vcf):
	index=0
	mvf_filter=""
	mvc_filter=""
	with open(vcf, "r") as f:
		for line in f:
			if "#" == line[0]:
				print(line.rstrip())
			if "##FILTER=<ID=MVF" in line:
				mvf_filter = line.split(",")[0].split("=")[2]
			if "##FILTER=<ID=MVC" in line:
				mvc_filter = line.split(",")[0].split("=")[2]
	f.close()
	return  mvf_filter, mvc_filter

def getSampleIndex(vcf, sampleName):
	index = 0
	with open(vcf, "r") as f:
		for line in f:
			if "#CHROM" in line:
				fields = line.rstrip().split("\t")
				index = fields.index(sampleName)
	f.close()
	return index	

def updateFilter(line, tumorIndex, normalIndex, keep, minvaf, mvf_filter_name, min_alt, mvc_filter_name):
	try:
		fields = line.split("\t")
		formatString = fields[8]
		adIndex = formatString.split(":").index("AD")
		afIndex = formatString.split(":").index("AF")
		tumorString = fields[tumorIndex]
		normalString = fields[normalIndex]
		alt_ad = int(tumorString.split(":")[adIndex].split(",")[1])
		af = float(tumorString.split(":")[afIndex])
		normal_af = float(normalString.split(":")[afIndex])
		if af < float(minvaf): 
			if "PASS" == fields[6]:
				fields[6] = mvf_filter_name
			elif mvf_filter_name not in fields[6]:
				fields[6] = fields[6] + ";" + mvf_filter_name
		if alt_ad < int(min_alt):
			if "PASS" == fields[6]:
				fields[6] = mvc_filter_name
			elif mvc_filter_name not in fields[6]:
				fields[6] = fields[6] + ";" + mvc_filter_name
		
		if keep and "whitelist" in line: # All whitelist calls are PASSed if keep parameter is set, regardless of past filtering (except for gnomAD filter)
			if 'filter_vep_fail' not in fields[6]:
				fields[6] = "PASS"
		# Do not pass any calls that have AF > 0.1%. Lower AF is assumed to be an error (Less than 1 read in 1000)
		if normal_af > 0.001:
			if fields[6] == "PASS":
				fields[6] = "germline"
			else:
				fields[6] = fields[6] + ";germline"
		updatedLine = "\t".join(fields)
		return updatedLine	
	except:
		return line
	

def applyFilter(vcf, samples, keep, minvaf, minalt):
	tumorName = samples.split(",")[0]
	normalName = samples.split(",")[1]
	if tumorName == normalName: # Unable to correctly distinguish columns
		printOriginal(vcf)
	else:
		mvf_filter_name, mvc_filter_name = printHeader(vcf)
		tumorIndex = getSampleIndex(vcf, tumorName)
		normalIndex = getSampleIndex(vcf, normalName)
		with open(vcf, "r") as f:
			for line in f:
				if "#" != line[0]:
					newline = updateFilter(line.rstrip(), tumorIndex, normalIndex, keep, minvaf, mvf_filter_name, minalt, mvc_filter_name)
					print(newline)	
		f.close()
		
def main(args):
	parser = argparse.ArgumentParser(description="Apply filters")
	parser.add_argument("-v", dest="vcf", action="store", help="Input vcf")
	parser.add_argument("-s", dest="samples", action="store", help="Comma-separated list in format: TumorSampleName,NormalSampleName")
	parser.add_argument("--filter", default=False, action="store_true", help="Apply filter. If flag is missing, output vcf is same as input vcf")
	parser.add_argument("--keep", default=False, action="store_true", help="Set all whitelist calls to PASS, except those that failed gnomAD filter")
	parser.add_argument("--mvf", dest="minvaf", action="store", help="Minimum vaf")
	parser.add_argument("--minalt", dest="minalt", action="store", help="Minimum alt support")
	args = parser.parse_args()
	
	if args.filter:
		applyFilter(args.vcf, args.samples, args.keep, args.minvaf, args.minalt)
	else:
		printOriginal(args.vcf)

if __name__ == '__main__':
	sys.exit(main(sys.argv))
