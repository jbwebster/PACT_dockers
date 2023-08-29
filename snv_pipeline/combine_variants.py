#!/usr/bin/python

'''
 Modified/custom implementation based on GATK3.6's CombineVariants
 Combines all SNVs/Indels from provided vcfs into a single file
 Assumes each input vcf comes from a different tool, with file name = .*<tool>_filtered.* , as is the case
 in the CWL workflow this is used for
'''

import argparse
import gzip
import re
import sys

class Call:
	def __init__(self, chrm, pos, ref, alt, f, tool, a, b):
		self.chrm = chrm
		self.pos = int(pos)
		self.ref = ref
		self.alt = alt
		self.filter = f.split(";")
		self.tool = tool
		self.a = a # sample a
		self.b = b # sample b

	def __eq__(self, other):
		return ((self.chrm, self.pos, self.ref, self.alt) == (other.chrm, other.pos, other.ref, other.alt))

	def add_tool(self, new_tool, f):
		if f != "PASS":
			new_tool = "filteredin" + new_tool
		self.tool = self.tool + "," + new_tool

	def add_filter(self, f):
		f = f.split(";")
		if "PASS" in f:
			self.filter = ["PASS"]
		elif "PASS" in self.filter:
			return
		else:
			for x in f:
				if x not in self.filter:
					self.filter.append(x)

	def to_string(self):
		filter_string = ";".join(self.filter)
		out = [self.chrm, str(self.pos), ".", str(self.ref), str(self.alt), ".", filter_string, "set=" + self.tool, "GT", str(self.a), str(self.b)]
		return "\t".join(out)

def define_parser():
	parser = argparse.ArgumentParser('combine-variants')
	parser.add_argument('-i', dest="input",  action='append')
	parser.add_argument('-o', dest="output", action="store")
	return parser

def update_headers(headers, line):
	if "##fileformat=" in line:
		return headers
	elif "##FILTER=" in line:
		if line not in headers:
			headers.append(line)
		return headers
	elif "##FORMAT=<ID=GT" in line:
		if line not in headers:
			headers.append(line)
		return headers
	elif "##INFO" in line:
		return headers
	elif "##contig" in line:
		if line not in headers:
			headers.append(line)
		return headers
	return headers

def parse_gt(form, sample):
	formats = form.split(":")
	sample_info = sample.split(":")
	index = formats.index("GT")
	return sample_info[index]
	

def parse_call(line):
	fields = line.split("\t")
	chrm = fields[0]
	pos = fields[1]
	ref = fields[3]
	alt = fields[4]
	f = fields[6]
	a = parse_gt(fields[8], fields[9])
	b = parse_gt(fields[8], fields[10])
	return chrm, pos, ref, alt, f, a, b
	

def update_calls(calls, line, tool):
	chrm, pos, ref, alt, f, a, b = parse_call(line)
	call = Call(chrm, pos, ref, alt, f, tool, a, b)
	found = False
	# Not a very good way of doing this..
	for i in range(0, len(calls)):
		prev_call = calls[i]
		if call == prev_call:
			prev_call.add_tool(tool, f)
			prev_call.add_filter(f)
			calls[i] = prev_call
			found = True
	if not found:
		calls.append(call)	
	return calls
		

def parse_path(path):
	pattern = 'inputs/.*/(.*)_filtered'
	return re.findall(pattern, path)[0]

def read_inputs(paths):
	# Assumed to be 5 inputs
	# (strelka, pindel, varscan, mutect, and whitelist)
	headers = ["##fileformat=VCFv4.2"]
	calls = []
	final_header = ""
	for path in paths:
		tool = parse_path(path)
		with gzip.open(path, 'rt') as f:
			for line in f:
				line = line.strip()
				if "##" in line:
					headers = update_headers(headers, line)
				elif "#CHROM" in line:
					final_header = line
				else:
					calls = update_calls(calls, line, tool)
		f.close()
	headers.append('##INFO=<ID=set,Number=.,Type=String,Description="Original source of call. May include variant caller(s) or whitelist as a source">')
	headers.append(final_header)
	return headers, calls	

def write_output(headers, calls, outfile):
	calls = sorted(calls, key=lambda x: (x.chrm, x.pos))
	with open(outfile, "w") as f:
		for header in headers:
			s = header + "\n"
			f.write(s)
		for call in calls:
			s = call.to_string() + "\n"
			f.write(s)
	f.close()

def main(args_input = sys.argv[1:]):
	parser = define_parser()
	args = parser.parse_args(args_input)
	input_paths = args.input
	headers, calls = read_inputs(input_paths)
	write_output(headers, calls, args.output)	
	


if __name__ == '__main__':
	sys.exit(main())



