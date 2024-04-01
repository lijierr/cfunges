#####################################################
#													#
# heatmap.py - plot identified genes as heatmap  	#
#													#
#####################################################

import matplotlib
matplotlib.use('Agg')
import matplitlib.pyplot as plt
import seaborn as sns

class heatmap:

	shape = (10, 6)
	print(shape)

	def __init__(self, abd, genes, annot_symbol=None):
		self.abd = abd
		self.genes = genes
		self.annot_symbol = annot_symbol
		
	def plot(self, outfig):
		self.abd = pd.read_csv(self.abd, sep='\t', header=0, index_col=0)
		self.genes = pd.read_csv(self.genes, sep='\t', header=0, index_col=0)
		plt.subplot2grid(shape, (0,0), rowspan=self.abd.shape[0], colspan=self.abd.shape[1])
		ax = sns.heatmap(self.abd, annot=True, cmap='RdBu_r',
						linecolor='k', linewidth=2, 
						xticklables=True, cbar=False)

		ax.xaxis.tick_top()
		plt.xticks(rotation=90)
		
		plt.subplot2grid(shape, (0,1), rowspan=self.genes.shape[0], colspan=self.genes.shape[1])
		if self.annot_symbol:
			if self.annot_symbol == 'numbers':
				ax = self._plot_symbol_number()
			else:
				ax = self._plot_symbol_symbol()
		else:
			ax = self._plot_no_annot()
		ax.xaxis.tick_top()
		ax.set_ylabel(' ') # to make sure no y axis label showed up
		plt.xticks(rotation=90)
		plt.savefig(outfig, dpt=800, bbox_inches='tight', pad_inches=.1)
	
	def _plot_symbol_number(self):
		'''plot here is you want to annotate heatmap with numbers.'''
		self.annot = self.genes.apply(lambda x: pd.Series(x)
								.apply(lambda y: 'NaN' if not y else y))
		ax = sns.heatmap(self.genes, annot=self.annot, cmap='RdBu',
						linecolor='k', linewidth=1.5,
						xticklabels=True, yticklabels=False,
						cbar=False)
		return ax

	def _plot_symbol_symbol(self):
		'''plot here if you want to annotate heatmap with symbol.'''
		self.annot = self.genes.apply(lambda x: pd.Series(x)
								.apply(lambda y: 'NaN' if not y else self.annot_symbol))
		ax = sns.heatmap(self.genes, annot=self.annot, fmt='s', cmap='RdBu',
						linecolor='k', linewidth=1.5,
						xticklabels=True, yticklabels=False,
						cbar=False)
		return ax

	def _plot_no_annot(self):
		'''plot here is you do not want to annotate heatmap.'''
		ax = sns.heatmap(self.genes, annot=False, cmap='RdBu',
						linecolor='k', linewidth=1.5,
						xticklabels=True, yticklabels=False,
						cbar=False)
		return ax


#	deprecated
#	self.symbol_annot and self._execute_heatplot function because too many if else of self.annot_symbol
#	def symbol_annot(self):
#		if self.annot_symbol == 'numbers':
#			self.annot = self.genes.apply(lambda x: pd.Series(x)
#									.apply(lambda y: 'NaN' if not y else y))
#		else:
#			self.annot = self.genes.apply(lambda x:pd.Series(x)
#									.apply(lambda y: 'NaN' if not y else annot_symbol))

#	def _execute_heatplot(self):
#		if self.annot_symbol == 'numbers':
#			ax = sns.heatmap(self.genes)
	
		

