#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'CC BY-NC-SA 4.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import os
import sys

from biolib.seq_io import read_fasta 
from GPMsDB_dbtk.common import (checkFileExists)

class Mw(object):
  def __init__(self):
      self.out_file_name = 'mw.txt'
      self.out_file_name2 = 'mw_s.txt'

  def run(self, aaFile):
      aawa = {
	  	'A' : 71.0788,  # alanine
	  	'R' : 156.1875, # arginine
	  	'D' : 115.0886, # aspartic acid
	  	'N' : 114.1038, # asparagine
	  	'C' : 103.1388, # cysteine
	  	'E' : 129.1155, # glutamic acid
	  	'Q' : 128.1307, # glutamine
	  	'G' : 57.0519,  # glycine
	  	'H' : 137.1411, # histidine
	  	'I' : 113.1594, # isoleucine
	  	'L' : 113.1594, # leucine
	  	'K' : 128.1741, # lysine
	  	'M' : 131.1926, # methionine
	  	'F' : 147.1766, # phenylalanine
	  	'P' : 97.1167,  # proline
	  	'S' : 87.0782,  # serine
	  	'T' : 101.1051, # threonine
	  	'W' : 186.2132, # tryptophan
	  	'Y' : 163.1760, # tyrosine
	  	'V' : 99.1326,  # valine
	      'U' : 150.0388, # selenocysteine
	      'O' : 237.3018 # pyrrolysine
      }
      # a list of molecular weights (monoisotopic) of the basic amino acids
      aawm = {
	  	'A' : 71.03711,  # alanine
	  	'R' : 156.10111, # arginine
	  	'D' : 115.02694, # aspartic acid
	  	'N' : 114.04293, # asparagine
	  	'C' : 103.00919, # cysteine
	  	'E' : 129.04259, # glutamic acid
	  	'Q' : 128.05858, # glutamine
	  	'G' : 57.02146,  # glycine
	  	'H' : 137.05891, # histidine
	  	'I' : 113.08406, # isoleucine
	  	'L' : 113.08406, # leucine
	  	'K' : 128.09496, # lysine
	  	'M' : 131.04049, # methionine
	  	'F' : 147.06841, # phenylalanine
	  	'P' : 97.05276,  # proline
	  	'S' : 87.03203,  # serine
	  	'T' : 101.04768, # threonine
	  	'W' : 186.07931, # tryptophan
	  	'Y' : 163.06333, # tyrosine
	  	'V' : 99.06841,  # valine
	      'U' : 150.953636, # selenocysteine
	      'O' : 237.147727 # pyrrolysine
      }
      checkFileExists(aaFile)

      file_dir, filename = os.path.split(aaFile)

      output_file = os.path.join(file_dir, self.out_file_name)

      fasta_sequences = read_fasta(aaFile)

      fout = open(output_file, 'w')
      fout.write('#Gene Id\taverage MH+\tmonoisotopic MH+\tSequence\n')
      ms_dic = {}

      for k in fasta_sequences.keys():
      	if fasta_sequences[k][0] == 'M':
            if fasta_sequences[k][1] == 'G' or fasta_sequences[k][1] == 'A' or fasta_sequences[k][1] == 'S' or fasta_sequences[k][1] == 'P' or fasta_sequences[k][1] == 'V' or fasta_sequences[k][1] == 'T' or fasta_sequences[k][1] == 'C':
              wa = 18.01056+1.00728-131.1926
              wm = 18.01056+1.00728-131.04049
            else:
              wa = 18.01056+1.00728
              wm = 18.01056+1.00728
      	else:
              wa = 18.01056+1.00728
              wm = 18.01056+1.00728

      	for aa in fasta_sequences[k]:
              if aa in aawa:
              	wa += aawa[aa]
              if aa in aawm:
              	wm += aawm[aa]

      	fout.write('%s\t%s\t%s\t%s\n' % (k, wa, wm, fasta_sequences[k]))
      	ms_dic[k] = wa

      fout.close

      result = sorted(ms_dic.items(), key=lambda x:x[1], reverse=False)
      key = dict(result)
      output_file2 = os.path.join(file_dir, self.out_file_name2)
      fout2 = open(output_file2, 'w')
      for j in key:
      	if 2000 < ms_dic[j] < 15000:
      		fout2.write('%s\t%s\n' % (j, ms_dic[j]))
      	else:
      	 	continue

      return ms_dic
