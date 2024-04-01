
#########################################################################
#																		#
# diamond.py - map querys to references, and give qulified results back.#
#																		#
#########################################################################

from .system_utils import check_software,file_utils
import subprocess
import pandas

class diamond:
	def __init__(self, query, subject, outdir):
		check_software('diamond')
		self.query = query
		self.subject = subject
		file_utils.check_file_exist([self.query, self.subject])
		self.outfile = outfile
		self.targs = targs
	
	def filter_diamond(self, targs):
		

	def _make_sure_diamond_db(self):
		'''Make diamond db'''
		if not os.path.isfile(self.subject + '.dmnd'):
			cmd = ['diamond', 'makedb', '--in', self.subject, '-d', self.subject]
			proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			proc.wait()
	
	def _diamond_map(self):
		'''Pileup diamond mapping'''
		_prefix = os.path.splitext(os.path.basename(self.query))[0]
		_suffix = os.path.splittext(os.path.basename(self.subject))[0]
		self.diamond_out = os.path.join(self.outdir, _prefix+'_diamond_to_'+_suffix+'.6')
		cmd = ['diamond', 'blastp', '-q', self.query, '-d', self.subject,
			'--outfmt', 6, '--strand', 'both', '--max-target-seqs', 10,
			'-o', self.diamond_out]
		proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		proc.wait()

	def pileup(self):
		if not os.path.isfile(self.subject + '.dmnd'):
			self._make_sure_diamond_db()
		self._diamond_map()

class diamond_parser:
	'''Parse diamond result (format 6)'''
	def __init__(self, m6, targs):
		self.m6 = m6
		self.targs = targs
		self.parsed_diamond = pandas.read_csv(self.targs.diamond_out, sep='\t', header=None, index_col=0)
		self.parsed_diamond.columns = ['sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
										'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'] ## add more columns to add length
		if self.targs.identity:self._identity()
		if self.targs.evalue:self._evalue()
		if self.targs.aln_length:self._aln_length()

	def _get_seq_length(self):
		from Bio import SeqIo
		seqs = SeqIO.to_dict(SeqIO.parse(self.query, 'fasta'))
		for i in seqs:
			seqs[i] = len(seqs[i])
		return pd.DataFrame.from_dict(seqs, orient='index', columns=['length'])

	def _qaln_ratio(self):
		'''drop alignments that not covered enough ratio of query length.'''
		qlens = self._get_seq_length()
		qaln_length = abs(self.parsed_diamond.qend-self.parsed_diamond.qstart+1)/qlens.loc[self.parsed_diamond.index].length

	def _aln_length(self):
		'''drop alignments that align length not meet cutoff.'''
		self.parsed_diamond = self.parsed_diamond[self.parsed_diamond.length>=self.targs.length]

	def _bitscore(self):
		'''drop alignments that not meet the bitscore cutoff.'''
		self.parsed_diamond = self.parsed_diamond[self.parsed_diamond.bitscore>=self.targs.bitscore]

	def _evalue(self):
		'''drop alignments that not meet evalue cutoff.'''
		self.parsed_diamond = self.parsed_diamond[self.parsed_diamond.evalue<=self.targs.evalue]

	def _identity(self):
		'''drop alignments that not meet identity cutoff.'''
		self.parsed_diamond = self.parsed_diamond[self.parsed_diamond.pident>=self.targs.identity]


