#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'CC BY-NC-SA 4.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import logging
import pickle

from GPMsDB_dbtk.common import checkFileExists
from GPMsDB_dbtk.defaultValues import DefaultValues


class Db(object):
  def __init__(self):
      self.db_file_r = DefaultValues.CUSTOM_LIST_R
      self.db_file_o = DefaultValues.CUSTOM_LIST_O
      self.db_file_genes = DefaultValues.CUSTOM_LIST_GENES
      self.db_file_names = DefaultValues.CUSTOM_LIST_NAME
      self.db_file_tax = DefaultValues.CUSTOM_LIST_TAX
      checkFileExists(self.db_file_r)
      checkFileExists(self.db_file_o)
      checkFileExists(self.db_file_genes)
      checkFileExists(self.db_file_names)
      checkFileExists(self.db_file_tax)
      self.logger = logging.getLogger('GPMsDB_tk')
      self.ribo_db = {}
      self.others_db = {}
      self.genes_db = {}
      self.tax_db = {}

  def list(self):
      with open(self.db_file_r, 'rb') as f:
         try:
             ribosomals = pickle.load(f)
             self.logger.info('[db_list] ' + str(len(ribosomals.keys())) + " entries found in the custom db")
             for i in ribosomals.keys():
                 if i == "":
                     pass
                 else:
                     self.logger.info(i)
         except EOFError:
             self.logger.info("no entry found in the custom db")

  def add(self, peakFile):
      checkFileExists(peakFile)
      Db.checkDb(self)

      genes = {}
      names = {}
      tax = {}
      ribosomals = {}
      others = {}

      a = 0
      for line in open(peakFile):
          if line.rstrip() == "":
              break
          if "Genome Id" in line:
              continue

          a += 1
          element = line.split("\t")
          id = "GCC_" + element[0].rstrip()
          ribosomals[id] = []
          others[id] = []
          if len(element) >= 7:
              genes[id] = int(element[1]) + int(element[2])
              names[id] = element[5]
              tax[id] = element[6].rstrip()
          elif len(element) == 6:
              genes[id] = int(element[1]) + int(element[2])
              names[id] = element[5].rstrip()
          else:
              genes[id] = int(element[1]) + int(element[2])

          item_r = element[3].replace('"',"").split(",")
          for n in item_r:
              ribosomals[id].append(n)
          item_o = element[4].replace('"',"").split(",")
          for j in item_o:
              others[id].append(j)

      Db.loadDb(self)

      for i in ribosomals.keys():
          if i in self.ribo_db.keys():
              self.logger.info('The same genome id was found in the new list')
              self.logger.info('Changing the id is reccomended for the genome: ' + i)
              break
      for j in ribosomals.keys():
          self.ribo_db[j] = ribosomals[j]
          self.others_db[j] = others[j]
          self.genes_db[j] = genes[j]
          try:
              self.names_db[j] = names[j]
          except:
              self.names_db[j] = ""
          try:
              self.tax_db[j] = tax[j]
          except:
              self.tax_db[j] = ""

      Db.dumpDb(self)

      self.logger.info(str(a) + " entries found and added in the custom db")


  def remove(self, accessions):
      Db.checkDb(self)
      Db.loadDb(self)

      list = []
      acces = accessions.split(",")
      for name in acces:
          name_fine = name.rstrip().lstrip()
          if name_fine in self.ribo_db.keys():
              list.append(name_fine)
          else:
              self.logger.info('Id is not found in the custom database: ' + name_fine)

      self.logger.info(str(len(list)) + ' genomes to be removed')

      for l in list:
          self.ribo_db.pop(l)
          self.others_db.pop(l)
          self.genes_db.pop(l)
          self.names_db.pop(l)
          self.tax_db.pop(l)

      Db.dumpDb(self)

  def checkDb(self):
      try:
          open(self.db_file_r, 'a')
          open(self.db_file_o, 'a')
          open(self.db_file_genes, 'a')
          open(self.db_file_names, 'a')
          open(self.db_file_tax, 'a')
      except IOError as e:
          self.logger.info("You do not seem to have permission to edit the custom db files")
          self.logger.info("Please try again with updated privileges. Error was:\n")
          self.logger.info(e)
          return False

  def loadDb(self):
      with open(self.db_file_r, 'rb') as f:
          self.ribo_db = pickle.load(f)
      with open(self.db_file_o, 'rb') as f:
          self.others_db = pickle.load(f)
      with open(self.db_file_genes, 'rb') as f:
          self.genes_db = pickle.load(f)
      with open(self.db_file_names, 'rb') as f:
          self.names_db = pickle.load(f)
      with open(self.db_file_tax, 'rb') as f:
          self.tax_db = pickle.load(f)

  def dumpDb(self):
      with open(self.db_file_r, mode='wb') as f:
          pickle.dump(self.ribo_db, f)
      with open(self.db_file_o, mode='wb') as f:
          pickle.dump(self.others_db, f)
      with open(self.db_file_genes, mode='wb') as f:
          pickle.dump(self.genes_db, f)
      with open(self.db_file_names, mode='wb') as f:
          pickle.dump(self.names_db, f)
      with open(self.db_file_tax, mode='wb') as f:
          pickle.dump(self.tax_db, f)
