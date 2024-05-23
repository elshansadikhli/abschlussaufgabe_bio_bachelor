#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import datei2 as dt2
import matplotlib.pyplot as plt
import numpy as np

class FASTAseq:
    with open('config.txt', 'r') as file:
        for line in file.readlines():
            if "fna" in line:
                path_fna = str(line.split('=')[1].rstrip())
            if "gff3" in line:
                path_gff3 = str(line.split('=')[1].rstrip())
    
    s1 = SeqIO.parse (path_fna,"fasta")
    for seq_rec in s1:
        full_sequence = str(seq_rec.seq)
        seqID = seq_rec.id
    gen1 = SeqRecord(Seq(full_sequence), seqID)
    
    
    def __init__(self):
        pass



    
    def creator_fasta (self):
        #Creating new Multifasta
        seq_wf_dict = {}
        if self.path_gff3 != None:
            gff3file = dt2.GFF3seq ()
            feat_dict = gff3file.creator_gff3 ()
            for k,v in feat_dict.items():
                seqrange = v.split (" ", 2)
                strand = str(seqrange[1])
                start = int(seqrange[0]) - 1
                end = int(seqrange[2])
                if strand == "+":
                    feat_seq = self.full_sequence[start:end]
                elif strand == "-":
                    feat_seq = self.full_sequence[start:end][::-1]
                else:
                    feat_seq = self.full_sequence[start:end]
                seq_wf_dict.update ({k : feat_seq})
            seqreclist = []
            aminolist = []
            for k,v in seq_wf_dict.items():
                seqrecs = SeqRecord(Seq(v), k)
                seqreclist.append(seqrecs)
                aminorecs = seqrecs.translate(id=True)
                aminolist.append(aminorecs)
            output_file_amino = "output_amino.fasta"
            output_file_nukl = "output_nukl.fasta"
            SeqIO.write(seqreclist, output_file_nukl, "fasta")
            SeqIO.write(aminolist, output_file_amino, "fasta")
        return seq_wf_dict
    
    def start_codone_usage(self):
        seq_wf_dict = FASTAseq.creator_fasta(self)
        list_of_startcods = []
        st_cod_analiser_dict = {}
        all_cod_alanizer_dict = {}
        for v in seq_wf_dict.values():
            startcodons = str(v)[:3]
            list_of_startcods.append(startcodons)
        st_cod_count_dict = dict(Counter(list_of_startcods))
        sum_st_cods = 0
        for val in st_cod_count_dict.values():
            sum_st_cods = sum_st_cods + int(val)
        for keys, vals in st_cod_count_dict.items():
            st_codshare = round((vals / sum_st_cods) * 100, 2)
            st_cod_analiser_dict.update({keys : str(str(vals)+" "+str(st_codshare)+"%")})
        all_cod_usage = [self.full_sequence[i:i + 3] for i in range(0, len(self.full_sequence), 3)]
        all_cod_count_dict = dict(Counter(all_cod_usage))
        sumallcods = 0
        for val in all_cod_count_dict.values():
            sumallcods = sumallcods + int(val)
        for keys, vals in all_cod_count_dict.items():
            all_codshare = round((vals / sumallcods)*100, 2)
            all_cod_alanizer_dict.update({keys : str(str(vals)+" "+str(all_codshare)+"%")})
        return st_cod_analiser_dict, all_cod_alanizer_dict
    
    def visualiser():
        plt.style.use('_mpl-gallery')

        # make data:
        x = 0.5 + np.arange(8)
        y = [4.8, 5.5, 3.5, 4.6, 6.5, 6.6, 2.6, 3.0]

        # plot
        fig, ax = plt.subplots()

        ax.bar(x, y, width=1, edgecolor="white", linewidth=0.7)

        ax.set(xlim=(0, 8), xticks=np.arange(1, 8),
            ylim=(0, 8), yticks=np.arange(1, 8))

        plt.show()

        























        


                
                


        

        
