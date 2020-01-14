
#########################################################################
#																		#
# diamond.py - map querys to references, and give qulified results back.#
#																		#
#########################################################################

from .system_utils import check_software,file_utils
import subprocess
import pandas

class diamond:
	def __init__(self, query, subject, outdir, **targs):
		check_software('diamond')
		self.query = query
		self.subject = subject
		file_utils.check_file_exist([self.query, self.subject])
		self.outfile = outfile

	def pileup(self):
		if not os.path.isfile(self.subject + '.dmnd'):
			self._mak_diamond_db()
		self._diamond_map()
		parse_diamond = diamond_parser(self.diamond_out)

	def _make_diamond_db(self):
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

class diamond_parser(self):
	'''Parse diamond result (format 6)'''
	def __init__(self, m6, **targs):
		self.m6 = m6
		self.targs = targs
		self.parsed_diamond = pandas.read_csv(self.targs.diamond_out, sep='\t', header=None, index_col=0)
		
	def _identity(self):
		self.parsed_diamond



