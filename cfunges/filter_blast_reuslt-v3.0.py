#!/usr/bin/env python
# File:filter_blast_reuslt-v3.0.py
# -*- coding: utf-8 -*-

"""
Contact: Jie LI (jie.li9197@gmail.com)
v1.0	201708	Choose blast result from m8(old blast) or m6(blast+)
v2.0	201809	Using pandas to be more faster
v3.0	201906	fixed the --top 1 not working bug, and more smarter
"""

import re
import sys
import argparse as ap
import gzip
import os
import subprocess as subp
import pandas as pd

def read_params(args):
	p = ap.ArgumentParser(description=__doc__)
	g = p.add_argument_group("Required arguments")
	ars = g.add_argument
	ars("--blast_m6", dest="blast_m6", required=True,
		help="m8 blast file input")
	ars("--outfile", dest="outfile", required=True,
		help="output file of chosen results")
	
	g = p.add_argument_group("Optional arguments")
	ars=g.add_argument
	ars("--uniq", dest="uniq", default=False, action="store_true",
		help="set to get a query with multi --top target but each targe are uniq")
	ars("--top", dest="top", default=None,
		help="choose the top num items, get all satisfied seqs if not set this")
	ars("--identity", dest="identity", default=40, type=int,
		help="identity threshold. [40]")
	ars("--align_length", dest="align_length", default=0, type=int,
		help="align length threshold.[0]")
	ars("--evalue", dest="evalue", default="1e-5",
		help="e-value threshold.[1e-5]")
	ars("--query_percent_length", dest="query_percent_length", default=False,
		help="min match percentage of length of query sequence, [0-100]")
	ars("--subject_percent_length", dest="subject_percent_length", default=False,
		help="min match percentage of length of subject sequence [0-100]")
	ars("--query", dest="query", default=False,
		help="query sequences file of blast if --query_percentage_length")
	ars("--subject", dest="subject", default=False,
		help="database sequences file of blast if --target_percentage_length")
	return vars(p.parse_args())

class filter_blast_result:
	def __init__(self, blast_m6, identity, evalue, align_length,
				 query_percent_length, subject_percent_length, 
				 query, subject, uniq, top, outfile):

		self.blast_m6 = blast_m6
		self.identity = float(identity)
		self.evalue = float(evalue)
		self.align_length = float(align_length)
		self.query_percent_length = query_percent_length
		self.subject_percent_length = subject_percent_length
		self.query = query
		self.subject = subject
		self.uniq = uniq
		self.top = top
		self.outfile = outfile
		self.filter_identity
		self.filter_evalue
		self.filter_identity()
		self.filter_evalue()
		self.check()	
		if self.align_length:self.filter_align_length()
		if self.query_percent_length:self.filter_query_percent_length()
		if self.subject_percent_length:self.filter_subject_percent_length()
		if self.uniq:self.filter_uniq()
		if self.top:self.filter_top()
		head = ['query', 'subject', 'identiy', 'length', 'mismatches', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
		if self.blast_m6.shape[0] >0 :
			self.blast_m6.to_csv(self.outfile, sep="\t", index=False, header=head)
		else:
			sys.exit("No blast result left! All just column names!\n")
		
	def check(self):
		if self.subject_percent_length and not self.subject:
			sys.exit(" Error: You set --subject_percent_length without providing --subject.")
		if self.query_percent_length and not self.query:
			sys.exit(" Error: You set --query_percent_length without providing --query.")

	def get_seq_length(self, fasta):
		from Bio import SeqIO
		sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
		for key in sequences:
			sequences[key] = len(sequences[key].seq)
		return sequences
	
	def filter_identity(self):
		print("Filtering identity not meet!")
		self.blast_m6 = self.blast_m6[self.blast_m6[2]>=self.identity]
		print("identity " + str(len(self.blast_m6)))

	def filter_evalue(self):
		print("Filtering evalue not meet!")
		self.blast_m6 = self.blast_m6[self.blast_m6[10]<=self.evalue]
		print("evalue " + str(len(self.blast_m6)))

	def filter_align_length(self):
		print("Filtering align lenght not meet!")
		self.blast_m6 = self.blast_m6[self.blast_m6[3]>=self.align_length]
		print("align_length " + str(len(self.blast_m6)))
	
	def filter_query_percent_length(self):
		print("Filtering algin length not meet the percentage of query cover!")
		self.query_percent_length = float(self.query_percent_length)
		self.query = self.get_seq_length(self.query)
		self.query_length = [self.query[i] for i in self.blast_m6[0]]
	#	self.cal_percent = abs(self.blast_m6[6]-self.blast_m6[7])*100/self.query_length
		self.cal_percent = abs(self.blast_m6[6]-self.blast_m6[7]+1)*100/self.query_length # 20191014, should +1, 这就是为什么当percent 为100的时候，全都判断为False的原因，不加1的时候，没有100%比对上的; 乌龙了，原来就是没有100%比对上的。。。
		#print(self.cal_percent)
		self.blast_m6 = self.blast_m6[self.cal_percent>=self.query_percent_length]   # 当percent为100的时候，这里的判断都是False， why? 20191014
		print("query_percent_length " + str(len(self.blast_m6)))

	def filter_subject_percent_length(self):
		print("Filtering align length not meet the percentage of subject cover!")
		self.subject_percent_length = float(self.subject_percent_length)
		self.subject = self.get_seq_length(self.subject)
		self.subject_length = [self.subject[i] for i in self.blast_m6[1]]
		self.cal_percent = abs(self.blast_m6[8]-self.blast_m6[9]+1)*100/self.subject_length
		self.blast_m6 = self.blast_m6[self.cal_percent>=self.subject_percent_length]
		print("subject_percent_length " + str(len(self.blast_m6)))

	def filter_uniq(self):
		uniq_queries = set(self.blast_m6[0])
		uniq_m6 = pd.DataFrame()
		for i in uniq_queries:
			i_blast_m6 = self.blast_m6[self.blast_m6[0]==i]
			for s in set(i_blast_m6[1]):
				uniq_m6 = uniq_m6.append(i_blast_m6[i_blast_m6[1]==s].iloc[[0]])
		self.blast_m6 = uniq_m6
		print("uniq " + str(len(self.blast_m6)))		

	def filter_top(self):
		self.top = int(self.top)
		top_m6 = pd.DataFrame()
		for i in set(self.blast_m6[0]):
			i_blast_m6 = self.blast_m6[self.blast_m6[0]==i]
			if len(i_blast_m6)>self.top:i_blast_m6 = i_blast_m6.iloc[0:self.top]
			top_m6 = top_m6.append(i_blast_m6)
		self.blast_m6 = top_m6	
		print("top " + str(len(self.blast_m6)))


if __name__ == "__main__":

	pars = read_params(sys.argv)
	if not os.path.getsize(pars['blast_m6']):sys.exit("Empty file input, please check!\n")
	blast_m6 = pd.read_csv(pars['blast_m6'], sep="\t", header=None, index_col=None)	
	print("Total lines " + str(len(blast_m6)))
	filter_blast_result(blast_m6, pars['identity'], pars['evalue'], pars['align_length'], pars['query_percent_length'], pars['subject_percent_length'], pars['query'], pars['subject'], pars['uniq'], pars['top'], pars['outfile'])


