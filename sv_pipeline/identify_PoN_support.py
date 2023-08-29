#!/usr/bin/python

import sys
import re
import argparse
import os.path

def removeNoise(vcf, percent):
	num_samples = 0
	with open(vcf, "r") as f:
		for line in f:
			line = line.rstrip()
			if "##" in line:
				print(line)
			elif "#CHROM" in line:
				print(line)
				cols =line.split("\t")
				samples = cols[9:]		
				num_samples = len(samples)
			else:
				fields = line.split("\t")	
				ad_index = fields[8].split(":").index("AD")
				dp_index = fields[8].split(":").index("DP")
				samples = fields[9:]
				supporting_samples = 0
				for s in samples:
					info = s.split(":")
					alt_reads = int(info[ad_index].split(",")[1]) * 1.0
					dp = int(info[dp_index]) * 1.0
					if alt_reads / dp > 0.001:
						supporting_samples += 1
				supporting_samples = supporting_samples * 1.0
				if (supporting_samples / num_samples) * 100 > percent:
						print(line)
	f.close()


def main(args):
	parser = argparse.ArgumentParser(description="Remove SNVs with support in a panel of normals")
	parser.add_argument("-v", "--vcf", action="store", help="Input VCF")
	parser.add_argument("-p", "--percent", action="store", default=100, type=int, help="Maximum percent of samples allows to support a variant. Default=10")
	args = parser.parse_args()

	removeNoise(args.vcf, args.percent)


if __name__ == '__main__':
	sys.exit(main(sys.argv))



