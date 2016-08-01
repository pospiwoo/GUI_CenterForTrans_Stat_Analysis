import os
import sys
import math
import re
import pickle
from sets import Set
import re
import scipy.stats
from scipy.stats import f_oneway
import numpy as np
import matplotlib.pyplot as plt
#import Heatmap
#import qvalue
#from statsmodels.sandbox.regression.predstd import wls_prediction_std
#import statsmodels
from statsmodels.sandbox.stats.multicomp import multipletests

#from rpy2.robjects.packages import importr
#from rpy2.robjects.vectors import FloatVector
#python Translational_Proteomics_Cecnter_Pipeline.py Params.txt test.txt



################################################################################
#########################     Model Class        ###############################
################################################################################
class TransPipelineModel():
	def __init__(self,
		search_tool ,
		column_peptide ,
		column_protein_id ,
		column_gene_symbol ,
		column_phospho ,
		Quantification_method ,
		Analysis_level ,
		Num_replicate ,
		Num_channels ,
		Phospho_probability_threshold_ ,
		Phospho_site_specificity ,
		search_tool_type ,
		result_file_names_and_channels ,
		group_1_list,
		group_2_list):

		self.search_tool = search_tool
		self.column_peptide = int(column_peptide)
		self.column_protein_id = int(column_protein_id)
		self.column_gene_symbol = int(column_gene_symbol)
		self.column_phospho = int(column_phospho)
		self.Quantification_method = Quantification_method
		self.Analysis_level = Analysis_level
		self.Num_replicate = int(Num_replicate)
		self.Num_channels = int(Num_channels)
		self.Phospho_probability_threshold_ = float(Phospho_probability_threshold_)
		self.Phospho_site_specificity = Phospho_site_specificity
		self.search_tool_type = search_tool_type
		self.result_file_names_and_channels = result_file_names_and_channels
		self.group_1_list = group_1_list.split(",")
		self.group_2_list = group_2_list.split(",")
		
		self.header_write_flag = -1
		self.peptide_or_gene_master_table = {}
		self.peptide_or_gene_master_table_zeros_removed = {}
		self.peptide_or_gene_master_table_median_represent = {}
		self.statistics_table_fold_change = {}
		self.statistics_table_TTest = {}

	def	ParseTable(self):
		if self.Analysis_level.upper() == 'PEPTIDE' or self.Analysis_level.upper() == 'PROTEIN' or self.Analysis_level.upper() == 'GENE':
			if self.Analysis_level.upper() == 'PEPTIDE':
				key_str_col = self.column_peptide
			elif self.Analysis_level.upper() == 'PROTEIN':
				key_str_col = self.column_protein_id
			elif self.Analysis_level.upper() == 'GENE':
				key_str_col = self.column_gene_symbol
			for i_rep in xrange(0,len(self.result_file_names_and_channels)):
				Result_File = open(self.result_file_names_and_channels[i_rep][0],'r')
				Result_raw = Result_File.readlines()
				count = 0
				if self.header_write_flag != -1:
					data = Result_raw[0].strip().replace('\'','').split('\t')
					for jjj in xrange(0,len(data)):
						print jjj, data[jjj]
				for j in range(1,len(Result_raw)):
					if Result_raw[j] == '':
						continue
					if Result_raw[j].find('DoNotUse')>0:
						continue
					if Result_raw[j].find('Not Found')>0:
						continue
					if Result_raw[j].find('NoQuanValues')>0:
						continue
					
					data = Result_raw[j].strip().replace('\"','').split('\t')
					
					#in case of gene name, we need to parse out gene name from protein accession
					if data[key_str_col].find('#') > -1:
						#XP_006710583.1#MACF1#23499#microtubule-actin cross-linking factor 1 isoform X5 [Homo sapiens]
						key_str = data[key_str_col].split('#')[1]
						#print key_str
					else:
						key_str = data[key_str_col]
					#print key_str, data[result_file_names_and_channels[i][1]], data[result_file_names_and_channels[i][2]], data[result_file_names_and_channels[i][3]]

					if self.peptide_or_gene_master_table.has_key(key_str):
						self.AddToExistingTableEntry(key_str, data, i_rep)
					else:
						self.AddNewTableEntry(key_str, data, i_rep)
				Result_File.close()

			#for jj in self.peptide_or_gene_master_table:
			#	print jj, self.peptide_or_gene_master_table[jj]





	def	AddToExistingTableEntry(self, key_str, data, i_rep):
		#print 'before', key_str, self.peptide_or_gene_master_table[key_str]
		for i in xrange(0,self.Num_channels):
			tmp_col = self.result_file_names_and_channels[i_rep][i+1]
			#print self.peptide_or_gene_master_table[key_str][i_rep][i], '+' , float(data[tmp_col])
			#print self.peptide_or_gene_master_table[key_str][i_rep][i], key_str, i_rep, i, data[tmp_col], tmp_col
			if tmp_col == -1:
				continue
			if data[tmp_col].strip() == '' or data[tmp_col].strip() == "":
				self.peptide_or_gene_master_table[key_str][i_rep][i] += 0.0
			else:
				self.peptide_or_gene_master_table[key_str][i_rep][i] += float(data[tmp_col])
		#print 'after ', key_str, self.peptide_or_gene_master_table[key_str]

	def	AddNewTableEntry(self, key_str, data, i_rep):
		######### peptide or gene master table ##########
		######### rep_1 channel 1 2 3 ##########
		######### rep_2 channel 1 2 3 ##########
		######### rep_3 channel 1 2 3 ##########
		tmp_replate_abun_table = []
		for i in xrange(0,self.Num_replicate):
			tmp_list = []
			for j in xrange(0,self.Num_channels):
				tmp_list.append(0.0)
			tmp_replate_abun_table.append(tmp_list)
		self.peptide_or_gene_master_table[key_str] = tmp_replate_abun_table
		
		for i in xrange(0,self.Num_channels):
			tmp_col = self.result_file_names_and_channels[i_rep][i+1]
			#print tmp_col, data[tmp_col]
			if tmp_col == -1:
				self.peptide_or_gene_master_table[key_str][i_rep][i] = 1.0
				continue
			if data[tmp_col] == '':
				self.peptide_or_gene_master_table[key_str][i_rep][i] = 0.0
			else:
				self.peptide_or_gene_master_table[key_str][i_rep][i] = float(data[tmp_col])




################################################################################
#######################     Controller Class        ############################
################################################################################
class TransPipelineController:
	def __init__(self, TPC_model_OBJ):
		TPC_model_OBJ.ParseTable()
		self.TPC_view = TransPipelineView()

	##############      Main processing function     ######################
	def process(self, TPC_model_OBJ):
		############# Remove zero entries #################
		print "before removing zero entries", len(TPC_model_OBJ.peptide_or_gene_master_table)
		self.TPC_view.PrintTable(TPC_model_OBJ.peptide_or_gene_master_table, '1_merged_table.txt')
		if TPC_model_OBJ.Quantification_method.upper() == "SILAC":
			self.RemoveNonIdentifiedEntriesSILAC(TPC_model_OBJ)
		elif TPC_model_OBJ.Quantification_method.upper() == "TMT":
			self.RemoveNonIdentifiedEntriesTMT(TPC_model_OBJ)
		print "after removing zero entries", len(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed)
		##PrintTableMedianLog2(peptide_or_gene_master_table_zeros_removed)
		self.TPC_view.PrintTable(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed, '2_zeros_removed.txt')
		
		
		############# Obtain median accross all replicates before normalization #################
		self.ObtainMedianAcrossReplicates(TPC_model_OBJ)		
		WFile = open('3_0_median_per_entry.txt','w')
		WFile.write('gene')
		for channel in xrange(0,TPC_model_OBJ.Num_channels):
			WFile.write('\tchannel_'+str(channel+1))
		WFile.write('\n')
		WFile.close()
		self.TPC_view.PrintTable(TPC_model_OBJ.peptide_or_gene_master_table_median_represent, '3_0_median_per_entry.txt')
		#for i in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table_median_represent)):
		#print scipy.stats.zscore(TPC_model_OBJ.peptide_or_gene_master_table_median_represent)

		############# Normanlization #################
		if TPC_model_OBJ.Quantification_method.upper() == "SILAC":
			self.SILACNormalization(TPC_model_OBJ)
		elif TPC_model_OBJ.Quantification_method.upper() == "TMT":
			self.TMTNormalization(TPC_model_OBJ)
		##PrintTable_light_med(peptide_or_gene_master_table_zeros_removed)
		self.TPC_view.PrintTable(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed, '3_normalized.txt')
		
		############# Obtain median accross all replicates #################
		self.ObtainMedianAcrossReplicates(TPC_model_OBJ)		
		WFile = open('3_1_normalized_and_median_per_entry.txt','w')
		WFile.write('gene')
		for channel in xrange(0,TPC_model_OBJ.Num_channels):
			WFile.write('\tchannel_'+str(channel+1))
		WFile.write('\n')
		WFile.close()
		self.TPC_view.PrintTable(TPC_model_OBJ.peptide_or_gene_master_table_median_represent, '3_1_normalized_and_median_per_entry.txt')
		
		############# Calculate Fold Change #################
		self.CalculateFoldChangeFromList(TPC_model_OBJ)
		self.TPC_view.PrintTableFoldChange(TPC_model_OBJ.statistics_table_fold_change, '5_fold_change.txt')
		
		############# Convert all abundance to log2 #################
		#self.Log2Transform(TPC_model_OBJ)
		#self.Log2TransformPlusOne(TPC_model_OBJ)
		#self.TPC_view.PrintTable(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed, '6_log2.txt')
		
		############# TTest #################
		#CalculateTTest(statistics_table_TTest,0,1)
		self.CalculateTTestFromList(TPC_model_OBJ)
		#TPC_model_OBJ.statistics_table_TTest.sort()
		self.TPC_view.PrintTablePValue(TPC_model_OBJ.statistics_table_TTest, '8_t_test.txt')



		############# Running R script for volcano plot #################
		RTableFile = open('9_R_table_volcano.txt','w')
		RTableFile.write('gene\tlog2FoldChange\tpvalue\n')
		for i in TPC_model_OBJ.statistics_table_TTest:
			#print i, statistics_table_fold_change[i]
			#RTableFile.write(i+"\t"+str(math.log(TPC_model_OBJ.statistics_table_fold_change[i],2))+'\t'+str(TPC_model_OBJ.statistics_table_TTest[i][1])+'\n')
			RTableFile.write(i+"\t"+str(TPC_model_OBJ.statistics_table_fold_change[i])+'\t'+str(TPC_model_OBJ.statistics_table_TTest[i][1])+'\n')
		RTableFile.close()
		#make R script
		R_txt = '''
		pdf(file="volcano_plot.pdf")
		res <- read.table(\"9_R_table_volcano.txt\", header=TRUE)
		'''
		# qvalue adjust
		R_txt += 'res$pvalue <- p.adjust(res$pvalue, "BH")\n'
		R_txt += '''
		with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10), ylim=c(0,8)))
		x  <- seq(-10, 10, 0.05)
		y <- x*0-log(0.05,10)
		lines(x,y,col="green")
		y <- x*0-log(0.01,10)
		lines(x,y,col="red")
		'''
		RScriptFile = open('9_R_script_volcano.r','w')
		RScriptFile.write(R_txt)
		RScriptFile.close()	
		cmd_line_str = "Rscript " + '9_R_script_volcano.r'
		os.system(cmd_line_str)



		############# Running R script for clustering and PCA #################
		#make R script for clustering
		R_txt = '''
		pdf(file="clustering_pca.pdf")
		library(pheatmap)
		res <- read.table(\"3_0_median_per_entry.txt\", header=TRUE)
		'''
		#res_sub <- res[,2:]
		R_txt += 'res_sub <- res[,2:' + str(TPC_model_OBJ.Num_channels+1) + ']\n'
		R_txt += 'res_sub <- log2(res_sub+1)\n'
		#pheatmap(res[,2:4])
		#R_txt += 'pheatmap(res[,2:' + str(TPC_model_OBJ.Num_channels+1) + '])\n'
		R_txt += 'pheatmap(res_sub)\n'
		#res_corr = cor(res)
		#R_txt += 'res_corr <- cor(res[,2:' + str(TPC_model_OBJ.Num_channels+1) + '])\n'
		R_txt += 'res_corr <- cor(res_sub)\n'
		R_txt += 'pheatmap(res_corr)\n'
		#make R script for PCA
		R_txt += 'library(ggbiplot)\n'
		R_txt += 'res_sub_t <- t(res_sub)\n'
		R_txt += 'res_pca <- prcomp(res_sub_t, center=TRUE, scale=TRUE)\n'
		R_txt += 'scores <- data.frame(rownames(res_sub_t), res_pca$x[,1:3])\n'
		R_txt += 'qplot(x=PC1, y=PC2, data=scores) + geom_text(label=rownames(res_sub_t))\n'
		RScriptFile = open('10_R_clustering_script.r','w')
		RScriptFile.write(R_txt)
		RScriptFile.close()	
		cmd_line_str = "Rscript " + '10_R_clustering_script.r'
		os.system(cmd_line_str)



	def	RemoveNonIdentifiedEntriesSILAC(self, TPC_model_OBJ):
		for i in TPC_model_OBJ.peptide_or_gene_master_table:
			cnt_zero = 0
			tmp_check_light_zero = -1
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				tmp_check_zero = sum(TPC_model_OBJ.peptide_or_gene_master_table[i][j])
				if tmp_check_zero == 0.0:
					cnt_zero += 1
				# if light channel abundence value is 0 we skip it, because we are going to normalize by light channel
				if float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][0]) == 0.0:
					tmp_check_light_zero = 1
			#if cnt_zero < Num_replicate-1 and tmp_check_light_zero == -1:
			if cnt_zero == 0 and tmp_check_light_zero == -1:
				TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed[i] = TPC_model_OBJ.peptide_or_gene_master_table[i]

	def	RemoveNonIdentifiedEntriesTMT(self, TPC_model_OBJ):
		for i in TPC_model_OBJ.peptide_or_gene_master_table:
			cnt_zero = 0
			tmp_check_light_zero = -1
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				tmp_check_zero = sum(TPC_model_OBJ.peptide_or_gene_master_table[i][j])
				if TPC_model_OBJ.result_file_names_and_channels[0][1] == -1 and tmp_check_zero == 1.0:
					cnt_zero += 1
				elif tmp_check_zero == 0.0:
					cnt_zero += 1
				# if the first channel abundence value is 0 we skip it, because we are going to normalize by the first channel
				if float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][0]) == 0.0:
					tmp_check_light_zero = 1
			#if cnt_zero < Num_replicate-1 and tmp_check_light_zero == -1:
			if cnt_zero == 0 and tmp_check_light_zero == -1:
				TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed[i] = TPC_model_OBJ.peptide_or_gene_master_table[i]

	def	SILACNormalization(self, TPC_model_OBJ):
		# Divide by light for each row
		for i in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				tmp_norm = float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][0])
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj] = float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj]) / tmp_norm

		# Divide by median for each col
		tmp_all_abun_per_channel = []
		for i in xrange(0,TPC_model_OBJ.Num_channels):
			tmp_empty = []
			tmp_all_abun_per_channel.append(tmp_empty)
		for i in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					tmp_all_abun_per_channel[jj].append(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj])
		median_per_channel = []
		for i in xrange(0,TPC_model_OBJ.Num_channels):
			tmp_med = np.median(tmp_all_abun_per_channel[i])
			#tmp_med = np.mean(tmp_all_abun_per_channel[i])
			median_per_channel.append(tmp_med)
		#print "median_per_channel",median_per_channel

		for i in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj] = float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj]) / float(median_per_channel[jj])

	def	TMTNormalization(self, TPC_model_OBJ):
		# Divide by the first channel for each row
		
		for i in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				tmp_norm = float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][0])
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj] = float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj]) / tmp_norm

		# Divide by median for each col
		tmp_all_abun_per_channel = []
		for i in xrange(0,TPC_model_OBJ.Num_channels):
			tmp_empty = []
			tmp_all_abun_per_channel.append(tmp_empty)
		for i in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					tmp_all_abun_per_channel[jj].append(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj])
		median_per_channel = []
		for i in xrange(0,TPC_model_OBJ.Num_channels):
			tmp_med = np.median(tmp_all_abun_per_channel[i])
			#tmp_med = np.mean(tmp_all_abun_per_channel[i])
			median_per_channel.append(tmp_med)
		#print "median_per_channel",median_per_channel

		for i in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj] = float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj]) / float(median_per_channel[jj])

		'''
		for i in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				#this is dividing by median across channel
				tmp_norm_list = []
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					tmp_norm_list.append(float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj]))
				tmp_norm = np.median(tmp_norm_list)

				if tmp_norm == 0.0:
					tmp_norm = 0.1
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj] = float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj]) / tmp_norm

		# Calculate median for each col
		tmp_all_abun_per_channel = []
		for i in xrange(0,TPC_model_OBJ.Num_channels):
			tmp_empty = []
			tmp_all_abun_per_channel.append(tmp_empty)
		for i in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					tmp_all_abun_per_channel[jj].append(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj])
		median_per_channel = []
		for i in xrange(0,TPC_model_OBJ.Num_channels):
			tmp_med = np.median(tmp_all_abun_per_channel[i])
			#tmp_med = np.mean(tmp_all_abun_per_channel[i])
			median_per_channel.append(tmp_med)

		for i in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for j in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i])):
				for jj in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[i][j])):
					TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj] = float(TPC_model_OBJ.peptide_or_gene_master_table[i][j][jj]) / float(median_per_channel[jj])
		'''




	def	Log2Transform(self, TPC_model_OBJ):
		for gene in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for rep in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[gene])):
				for channel in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[gene][rep])):
					TPC_model_OBJ.peptide_or_gene_master_table[gene][rep][channel] = math.log( float(TPC_model_OBJ.peptide_or_gene_master_table[gene][rep][channel]) , 2)

	def	Log2TransformPlusOne(self, TPC_model_OBJ):
		for gene in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for rep in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[gene])):
				for channel in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table[gene][rep])):
					TPC_model_OBJ.peptide_or_gene_master_table[gene][rep][channel] = math.log( float(TPC_model_OBJ.peptide_or_gene_master_table[gene][rep][channel])  + 1.0, 2)

	def	CalculateTTest(statistics_table_TTest, col_group_1, col_group_2, TPC_model_OBJ):
		######### peptide or gene master table ##########
		######### rep_1 channel 1 2 3 ##########
		######### rep_2 channel 1 2 3 ##########
		######### rep_3 channel 1 2 3 ##########
		for i in peptide_or_gene_master_table_zeros_removed:
			tmp_list_1 = []
			tmp_list_2 = []
			for j in xrange(0,len(peptide_or_gene_master_table_zeros_removed[i])):
				tmp_list_1.append(peptide_or_gene_master_table_zeros_removed[i][j][col_group_1])
				tmp_list_2.append(peptide_or_gene_master_table_zeros_removed[i][j][col_group_2])
			#print len(tmp_list_1), len(tmp_list_2), tmp_list_1, tmp_list_2
			ttest_val, two_tail_pvalue = scipy.stats.ttest_ind(tmp_list_1,tmp_list_2, equal_var = True)
			#ttest_val, two_tail_pvalue = scipy.stats.ttest_ind(tmp_list_1,tmp_list_2, equal_var = False)
			#print ttest_val, two_tail_pvalue, "~", tmp_list_1, tmp_list_2
			statistics_table_TTest[i] = [ttest_val, two_tail_pvalue]
			#statistics_table_TTest[i] = [ttest_val, -math.log(two_tail_pvalue,10)]


	def	CalculateTTestFromList(self, TPC_model_OBJ):
		######### peptide or gene master table ##########
		######### rep_1 channel 1 2 3 ##########
		######### rep_2 channel 1 2 3 ##########
		######### rep_3 channel 1 2 3 ##########
		for gene in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			tmp_list_1 = []
			tmp_list_2 = []
			for rep in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed[gene])):
				for channel in xrange(0,len(TPC_model_OBJ.group_1_list)):
					tmp_list_1.append(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed[gene][rep][int(TPC_model_OBJ.group_1_list[channel])-1])
				for channel in xrange(0,len(TPC_model_OBJ.group_2_list)):
					tmp_list_2.append(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed[gene][rep][int(TPC_model_OBJ.group_2_list[channel])-1])
			#print len(tmp_list_1), len(tmp_list_2), tmp_list_1, tmp_list_2
			ttest_val, two_tail_pvalue = scipy.stats.ttest_ind(tmp_list_1,tmp_list_2, equal_var = True)
			#ttest_val, two_tail_pvalue = scipy.stats.ttest_ind(tmp_list_1,tmp_list_2, equal_var = False)
			#print ttest_val, two_tail_pvalue, "~", tmp_list_1, tmp_list_2
			TPC_model_OBJ.statistics_table_TTest[gene] = [ttest_val, two_tail_pvalue]
			#statistics_table_TTest[i] = [ttest_val, -math.log(two_tail_pvalue,10)]

	def	CalculateFoldChange(statistics_table_fold_change,col_denominator,col_compare, TPC_model_OBJ):
		for gene in peptide_or_gene_master_table_zeros_removed:
			tmp_list_1 = []
			tmp_list_2 = []
			for rep in xrange(0,len(peptide_or_gene_master_table_zeros_removed[gene])):
				tmp_list_1.append(peptide_or_gene_master_table_zeros_removed[gene][rep][col_denominator])
				tmp_list_2.append(peptide_or_gene_master_table_zeros_removed[gene][rep][col_compare])

			tmp_median_1 = np.median(tmp_list_1)
			tmp_median_2 = np.median(tmp_list_2)
			#print tmp_list_1, tmp_median_1,  tmp_list_2, tmp_median_2,  tmp_median_2 / tmp_median_1
			statistics_table_fold_change[gene] = tmp_median_2 / tmp_median_1

	def	CalculateFoldChangeFromList(self, TPC_model_OBJ):
		for gene in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			tmp_list_1 = []
			tmp_list_2 = []
			for rep in xrange(0,len(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed[gene])):
				for group in xrange(0,len(TPC_model_OBJ.group_1_list)):
					tmp_list_1.append(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed[gene][rep][int(TPC_model_OBJ.group_1_list[group])-1])
				for group in xrange(0,len(TPC_model_OBJ.group_2_list)):
					tmp_list_2.append(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed[gene][rep][int(TPC_model_OBJ.group_2_list[group])-1])

			tmp_median_1 = np.median(tmp_list_1)
			tmp_median_2 = np.median(tmp_list_2)
			#print tmp_list_1, tmp_median_1,  tmp_list_2, tmp_median_2,  tmp_median_2 / tmp_median_1
			if tmp_median_1 == 0.0:
				tmp_median_1 = 0.1
			if tmp_median_2 == 0.0:
				tmp_median_2 = 0.1
			TPC_model_OBJ.statistics_table_fold_change[gene] = tmp_median_2 / tmp_median_1
			#if i == "EYA3":
			#	print "2222", tmp_list_1, tmp_list_2
			#	print "111111", tmp_median_1, tmp_median_2, statistics_table_fold_change[i]

	def	ObtainMedianAcrossReplicates(self, TPC_model_OBJ):
		######### peptide or gene median abundence table ##########
		######### gene channel 1 2 3 ##########
		for gene in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			tmp_list = []
			for channel in xrange(0,TPC_model_OBJ.Num_channels):
				tmp_list.append(0.0)
			TPC_model_OBJ.peptide_or_gene_master_table_median_represent[gene] = tmp_list
		for gene in TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed:
			for channel in xrange(0,TPC_model_OBJ.Num_channels):
				tmp_list = []			
				for rep in xrange(0,TPC_model_OBJ.Num_replicate):
					tmp_list.append(TPC_model_OBJ.peptide_or_gene_master_table_zeros_removed[gene][rep][channel])
				tmp_median = np.median(tmp_list)
				TPC_model_OBJ.peptide_or_gene_master_table_median_represent[gene][channel] = tmp_median


	def p_adjust_bh(p, TPC_model_OBJ):
	    '''Benjamini-Hochberg p-value correction for multiple hypothesis testing.'''
	    p = np.asfarray(p)
	    by_descend = p.argsort()[::-1]
	    by_orig = by_descend.argsort()
	    steps = float(len(p)) / np.arange(len(p), 0, -1)
	    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	    return q[by_orig]


	def calc_benjamini_hochberg_corrections(p_values, num_total_tests, TPC_model_OBJ):
	    '''
	    Calculates the Benjamini-Hochberg correction for multiple hypothesis
	    testing from a list of p-values *sorted in ascending order*.

	    See
	    http://en.wikipedia.org/wiki/False_discovery_rate#Independent_tests
	    for more detail on the theory behind the correction.

	    **NOTE:** This is a generator, not a function. It will yield values
	    until all calculations have completed.

	    :Parameters:
	    - `p_values`: a list or iterable of p-values sorted in ascending
	      order
	    - `num_total_tests`: the total number of tests (p-values)
	    '''
	    prev_bh_value = 0
	    for i, p_value in enumerate(p_values):
		bh_value = p_value * num_total_tests / (i + 1)
		# Sometimes this correction can give values greater than 1,
		# so we set those values at 1
		bh_value = min(bh_value, 1)

		# To preserve monotonicity in the values, we take the
		# maximum of the previous value or this one, so that we
		# don't yield a value less than the previous.
		bh_value = max(bh_value, prev_bh_value)
		prev_bh_value = bh_value
		yield bh_value


	def binary_search(array,key,imin,imax, TPC_model_OBJ):
	    if (imax < imin):
		return imax
	    else:
		imid = (imin+imax)/2
		if array[imid] > key:
		    return binary_search(array,key,imin,imid-1)
		elif array[imid] < key:
		    return binary_search(array,key,imid+1,imax)
		else:
		    return imid
		

	def removekey(d, key, TPC_model_OBJ):
	    r = dict(d)
	    del r[key]
	    return r


	def inserts(original, new, pos, TPC_model_OBJ):
		return original[:pos] + new + original[pos:]


	def CleanPeptideString(pep_str, TPC_model_OBJ):
		#print pep_str
		pep_str = pep_str.replace('0','')
		pep_str = pep_str.replace('1','')
		pep_str = pep_str.replace('2','')
		pep_str = pep_str.replace('3','')
		pep_str = pep_str.replace('4','')
		pep_str = pep_str.replace('5','')
		pep_str = pep_str.replace('6','')
		pep_str = pep_str.replace('7','')
		pep_str = pep_str.replace('8','')
		pep_str = pep_str.replace('9','')
		pep_str = pep_str.replace('+','')
		pep_str = pep_str.replace('.','')
		pep_str = pep_str.replace('?','_')
		pep_str = pep_str.replace('_','')
		pep_str = pep_str.replace('-','')
		pep_str = pep_str.replace('*','')
		#pep_str = inserts(pep_str,'.',1)
		#pep_str = inserts(pep_str,'.',-1)
		#print pep_str
		return pep_str







#####################################################################
########################### View class ##############################
#####################################################################
class TransPipelineView:
	def __init__(self):
		pass

	def	PrintTable(self, table_to_print, tmp_output_file_name):
		printTableFile = open(tmp_output_file_name,'a')
		#printTableFile.write('gene\tlog2FoldChange\tpvalue\n')
		for i in table_to_print:
			printTableFile.write(i)
			for j in xrange(0,len(table_to_print[i])):
				printTableFile.write("\t"+str(table_to_print[i][j]))
			printTableFile.write('\n')
		printTableFile.close()


	def	PrintTableFoldChange(self, table_to_print, tmp_output_file_name):
		printTableFile = open(tmp_output_file_name,'w')
		#printTableFile.write('gene\tlog2FoldChange\tpvalue\n')
		for i in table_to_print:
			printTableFile.write(i)
			printTableFile.write("\t"+str(table_to_print[i]))
			printTableFile.write('\n')
		printTableFile.close()



	def	PrintTablePValue(self, table_to_print, tmp_output_file_name):
		printTableFile = open(tmp_output_file_name,'w')
		printTableFile.write('gene\tpvalue\n')
		for i in table_to_print:
			printTableFile.write(i)
			printTableFile.write("\t"+str(table_to_print[i][1]))
			printTableFile.write('\n')
		printTableFile.close()


	def	PrintTableMedianLog2(table_to_print, TPC_model_OBJ):
		printTableFile = open('print_table_for_test.txt','w')
		printTableFile.write('gene\tLight\tMedium\tHeavy\n')
		for i in table_to_print:
			printTableFile.write(i)
			tmp_list_1 = []
			tmp_list_2 = []
			tmp_list_3 = []
			for j in xrange(0,len(table_to_print[i])):
				tmp_list_1.append(table_to_print[i][j][0])
				tmp_list_2.append(table_to_print[i][j][1])
				tmp_list_3.append(table_to_print[i][j][2])
			tmp_median_1 = np.median(tmp_list_1)
			tmp_median_2 = np.median(tmp_list_2)
			tmp_median_3 = np.median(tmp_list_3)

			printTableFile.write("\t"+str(math.log(tmp_median_1,2)))
			printTableFile.write("\t"+str(math.log(tmp_median_2,2)))
			printTableFile.write("\t"+str(math.log(tmp_median_3,2)))
			printTableFile.write('\n')
		printTableFile.close()


	def	PrintTable_light_med(table_to_print, TPC_model_OBJ):
		printTableFile = open('print_table_for_test.txt','w')
		#printTableFile.write('gene\tlog2FoldChange\tpvalue\n')
		for i in table_to_print:
			printTableFile.write(i)
			printTableFile.write("\t"+str(math.log(table_to_print[i][0][0],2)))
			printTableFile.write("\t"+str(math.log(table_to_print[i][1][0],2)))
			printTableFile.write("\t"+str(math.log(table_to_print[i][2][0],2)))
			printTableFile.write("\t"+str(math.log(table_to_print[i][0][1],2)))
			printTableFile.write("\t"+str(math.log(table_to_print[i][1][1],2)))
			printTableFile.write("\t"+str(math.log(table_to_print[i][2][1],2)))
			printTableFile.write('\n')
		printTableFile.close()









