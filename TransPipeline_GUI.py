#!/usr/bin/python
import os
import sys
import math
import re
import pickle
from sets import Set
#import scipy.stats
#from scipy.stats import f_oneway
#import numpy as np
#import matplotlib.pyplot as plt
#import Heatmap
#from statsmodels.sandbox.regression.predstd import wls_prediction_std
#import statsmodels
#from statsmodels.sandbox.stats.multicomp import multipletests
import tkMessageBox
import Tkinter, Tkconstants, tkFileDialog
import TransPipelineModules

class TransPipelineGUI(Tkinter.Frame):

	def __init__(self, root):

		self.input_pickle_file_name = ''
		self.header_write_flag = -1
		
		self.tmp_rep_file_name = ''
		self.params_file = ''
		self.output_file_name = 'out.txt'
		self.group_1_list = '1'
		self.group_2_list = '2'

		self.search_tool = ''
		self.column_peptide = 0
		self.column_protein_id = 0
		self.column_gene_symbol = 0
		self.column_phospho = 0
		self.Quantification_method = ''
		self.Analysis_level = ''
		self.Num_replicate = 0
		self.Num_channels = 0
		self.Phospho_probability_threshold_ = 75
		self.Phospho_site_specificity = ''

		self.search_tool_type = ''
		self.result_file_names_and_channels = []
		self.tmp_channels = []
		self.tmp_result_file_names = []
		#self.groups = []

		self.peptide_or_gene_master_table = {}
		self.peptide_or_gene_master_table_zeros_removed = {}

		Tkinter.Frame.__init__(self, root)

		# options for buttons
		#button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}

		# define options for opening or saving a file
		self.file_opt = options = {}
		#options['defaultextension'] = '.txt'
		#options['filetypes'] = [('all files', '.*'), ('text files', '.txt')]
		#options['initialdir'] = 'C:\\Research\\Software\\TransCenterPipeline\\GUI'
		#options['initialdir'] = 'C:\\Users\\Sunghee\\'
		#options['initialfile'] = 'myfile.txt'
		#options['parent'] = root
		#options['title'] = 'Result file'


		# define buttons
		button_row_ind = 0
		Tkinter.Label(self, text='Search Tool').grid(row=button_row_ind,column=0)
		self.search_tool_var = Tkinter.StringVar()
		Tkinter.Radiobutton(self, text='PD', variable=self.search_tool_var, value='PD', command=self.assign_search_tool).grid(row=button_row_ind,column=1)
		Tkinter.Radiobutton(self, text='MaxQuant', variable=self.search_tool_var, value='MaxQuant', command=self.assign_search_tool).grid(row=button_row_ind,column=2)
		Tkinter.Radiobutton(self, text='MSGF', variable=self.search_tool_var, value='MSGF', command=self.assign_search_tool).grid(row=button_row_ind,column=3)

		button_row_ind += 1
		Tkinter.Label(self, text='Analysis level').grid(row=button_row_ind,column=0)
		self.analysis_level_var = Tkinter.StringVar()
		Tkinter.Radiobutton(self, text='Peptied', variable=self.analysis_level_var, value='PEPTIDE', command=self.assign_analysis_level).grid(row=button_row_ind,column=1)
		Tkinter.Radiobutton(self, text='Gene', variable=self.analysis_level_var, value='GENE', command=self.assign_analysis_level).grid(row=button_row_ind,column=2)
		Tkinter.Radiobutton(self, text='Protein', variable=self.analysis_level_var, value='PROTEIN', command=self.assign_analysis_level).grid(row=button_row_ind,column=3)
		Tkinter.Radiobutton(self, text='Phospho', variable=self.analysis_level_var, value='PHOSPHO', command=self.assign_analysis_level).grid(row=button_row_ind,column=4)

		button_row_ind += 1
		Tkinter.Label(self, text='Peptide sequence column(0-based)').grid(row=button_row_ind,column=0)
		self.column_peptide_var = Tkinter.Spinbox(self, from_=0, to=500, width=5)
		self.column_peptide_var.grid(row=button_row_ind,column=1)
		
		button_row_ind += 1
		Tkinter.Label(self, text='Quantification method').grid(row=button_row_ind,column=0)
		self.Quantification_method_var = Tkinter.StringVar()
		Tkinter.Radiobutton(self, text='SILAC', variable=self.Quantification_method_var, value='SILAC', command=self.assign_quant).grid(row=button_row_ind,column=1)
		Tkinter.Radiobutton(self, text='TMT', variable=self.Quantification_method_var, value='TMT', command=self.assign_quant).grid(row=button_row_ind,column=2)

		button_row_ind += 1
		Tkinter.Label(self, text='Number of replicates').grid(row=button_row_ind,column=0)
		self.Num_replicate_var = Tkinter.Spinbox(self, from_=1, to=1000, width=5)
		self.Num_replicate_var.grid(row=button_row_ind,column=1)

		button_row_ind += 1
		Tkinter.Label(self, text='Number of different abundence channels').grid(row=button_row_ind,column=0)
		self.Num_channels_var = Tkinter.Spinbox(self, from_=1, to=50, width=5)
		self.Num_channels_var.grid(row=button_row_ind,column=1)

		button_row_ind += 1
		Tkinter.Button(self, text ='Submit', command = self.submit).grid(row=button_row_ind,column=0)

	def assign_quant(self):
		self.Quantification_method = self.Quantification_method_var.get()

	def assign_analysis_level(self):
		self.Analysis_level = self.analysis_level_var.get()

	def assign_search_tool(self):
		self.search_tool = self.search_tool_var.get()

	def submit(self):
		self.column_peptide = int(self.column_peptide_var.get())
		self.Num_replicate = int(self.Num_replicate_var.get())
		self.Num_channels = int(self.Num_channels_var.get())
		
		#root.quit()
		root.destroy()
		top = Tkinter.Toplevel()
		button_row_ind = 0

		if self.Analysis_level == 'PEPTIDE':
			print 'Column of peptide sequence', self.column_peptide

		elif self.Analysis_level == 'PROTEIN':
			Tkinter.Label(top, text='Protein ID column(0-based)').grid(row=button_row_ind,column=0)
			self.column_protein_id_var = Tkinter.Spinbox(top, from_=0, to=500, width=5)
			self.column_protein_id_var.grid(row=button_row_ind,column=1)
			
		elif self.Analysis_level == 'GENE':
			Tkinter.Label(top, text='Gene ID(symbol) column(0-based)').grid(row=button_row_ind,column=0)
			self.column_gene_symbol_var = Tkinter.Spinbox(top, from_=0, to=500, width=5)
			self.column_gene_symbol_var.grid(row=button_row_ind,column=1)

		elif self.Analysis_level == 'PHOSPHO':
			Tkinter.Label(top, text='Phospho sites column(0-based)').grid(row=button_row_ind,column=0)
			self.column_phospho_var = Tkinter.Spinbox(top, from_=0, to=500, width=5)
			self.column_phospho_var.grid(row=button_row_ind,column=1)
			self.Phospho_probability_threshold_ = 75
	#		self.Phospho_site_specificity = ''							
	#		Phospho_site_specificity	UNIQUE_PEPTIDE_SITES

		else:
			print 'Error! unknown Analysis level:', self.Analysis_level
		

		#create TMT/SILAC table
		button_row_ind += 1
		
		for i in xrange(0,self.Num_replicate):
			self.tmp_result_file_names.append('')
			tmp_list = []
			for j in xrange(0,self.Num_channels):
				tmp_list.append(Tkinter.IntVar())
			self.tmp_channels.append(tmp_list)

		if self.Quantification_method == 'SILAC':
			Tkinter.Label(top, text='Selected result files').grid(row=button_row_ind,column=1)
			for j in xrange(0,self.Num_channels):
				Tkinter.Label(top, text='Channel '+str(j+1)).grid(row=button_row_ind,column=j+2)
			button_row_ind += 1
			for i in xrange(0,self.Num_replicate):
				Tkinter.Label(top, text='Rep '+str(i+1)+' abundence columns(0-based)').grid(row=button_row_ind,column=0)
				Tkinter.Button(top, text='Rep file '+str(i+1), command=self.askRepResultFilename(i)).grid(row=button_row_ind,column=1)
				Tkinter.Label(top, text=os.path.basename(self.tmp_result_file_names[i])).grid(row=button_row_ind,column=1)
				tmp_list.append(str(self.tmp_rep_file_name))
				for j in xrange(0,self.Num_channels):

					self.tmp_channels[i][j] = Tkinter.Spinbox(top, from_=0, to=500, width=5)
					self.tmp_channels[i][j].grid(row=button_row_ind,column=j+2)

				button_row_ind += 1
			print 'Quantification_method:', self.Quantification_method
		elif self.Quantification_method == 'TMT':
			print 'Quantification_method:', self.Quantification_method
		else:
			print 'Error! unknown Quantification_method:', self.Quantification_method


		#Group 1
		button_row_ind += 1
		Tkinter.Label(top, text='Statistical test group 1').grid(row=button_row_ind,column=0)
		self.group_1_list_var = []
		for j in xrange(0,self.Num_channels):
			self.group_1_list_var.append(Tkinter.IntVar())
			Tkinter.Checkbutton(top, text='channel'+str(j+1), variable=self.group_1_list_var[j]).grid(row=button_row_ind,column=j+2)

		#Group 2
		button_row_ind += 1
		Tkinter.Label(top, text='Statistical test group 2').grid(row=button_row_ind,column=0)
		self.group_2_list_var = []
		for j in xrange(0,self.Num_channels):
			self.group_2_list_var.append(Tkinter.IntVar())
			Tkinter.Checkbutton(top, text='channel'+str(j+1), variable=self.group_2_list_var[j]).grid(row=button_row_ind,column=j+2)

		#Run button
		button_row_ind += 1
		Tkinter.Button(top, text ='Run', command = self.run_main).grid(row=button_row_ind,column=0)


##############################################

	def run_main(self):
		if self.Analysis_level == 'PROTEIN':
			self.column_protein_id = self.column_protein_id_var.get()
			print 'column_protein_id_var', self.column_protein_id_var.get()
		elif self.Analysis_level == 'GENE':
			self.column_gene_symbol = self.column_gene_symbol_var.get()
			print 'column_gene_symbol_var', self.column_gene_symbol_var.get()
		elif self.Analysis_level == 'PHOSPHO':
			self.column_phospho = self.column_phospho_var.get()
			print 'column_phospho_var', self.column_phospho_var.get()
			
		for i in xrange(0,self.Num_replicate):
			tmp_list = []
			tmp_list.append(self.tmp_result_file_names[i])
			for j in xrange(0,self.Num_channels):
				tmp_list.append(-1)
			self.result_file_names_and_channels.append(tmp_list)

		for i in xrange(0,self.Num_replicate):
			for j in xrange(0,self.Num_channels):
				self.result_file_names_and_channels[i][j+1] = int(self.tmp_channels[i][j].get())

		self.group_1_list = ''
		for j in xrange(0,self.Num_channels):
			if self.group_1_list_var[j].get() == 1:
				self.group_1_list = self.group_1_list + str(j+1)+','
		self.group_1_list = self.group_1_list.strip(',')
		
		self.group_2_list = ''
		for j in xrange(0,self.Num_channels):
			if self.group_2_list_var[j].get() == 1:
				self.group_2_list = self.group_2_list + str(j+1)+','
		self.group_2_list = self.group_2_list.strip(',')

#		for i in xrange(0,self.Num_replicate):
#			tmp_str = ''
#			for j in xrange(0,self.Num_channels):
#				tmp_str += ' ' + str(self.result_file_names_and_channels[i][j+1])
#			tmp_str += ' ' + os.path.basename(self.result_file_names_and_channels[i][0])
#			print tmp_str
#		print 'group_1_list', self.group_1_list
#		print 'group_2_list', self.group_2_list


		#Run Main
		self.TransPipelineRun()



	def askopenfilename(self):
		self.params_file = tkFileDialog.askopenfilename(**self.file_opt)

	def askRepResultFilename(self, tmp_ind):
		self.tmp_result_file_names[tmp_ind] = tkFileDialog.askopenfilename(title="Select Result file for Replicate "+str(tmp_ind+1))
		print 'File name Replicate', tmp_ind, ':', self.tmp_result_file_names[tmp_ind]




	def	TransPipelineRun(self):
		input_pickle_file_name = ''
		header_write_flag = -1

		############# init model controller #################
		TPM_model = TransPipelineModules.TransPipelineModel(
			self.search_tool ,
			self.column_peptide ,
			self.column_protein_id ,
			self.column_gene_symbol ,
			self.column_phospho ,
			self.Quantification_method ,
			self.Analysis_level ,
			self.Num_replicate ,
			self.Num_channels ,
			self.Phospho_probability_threshold_ ,
			self.Phospho_site_specificity ,
			self.search_tool_type ,
			self.result_file_names_and_channels ,
			self.group_1_list ,
			self.group_2_list )
		TPM_controller = TransPipelineModules.TransPipelineController(TPM_model)

		############# Run module controller #################
		TPM_controller.process(TPM_model)


if __name__=='__main__':
  root = Tkinter.Tk()
  TransPipelineGUI(root).pack()
  root.mainloop()




