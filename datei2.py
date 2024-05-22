#!/usr/bin/env python3
import pandas as pd
from collections import Counter


class GFF3seq:
    #Reading the file path from the config file
    with open('config.txt', 'r') as file:
        for line in file.readlines():
            if "gff3" in line:
                path_gff3 = line.split('=')[1]
                path = str(path_gff3.rstrip())
    #Parsing of GFF3 File, turning into DataFrame
    #Initially separating if whitespace
    column_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    df = pd.read_csv(path, sep='\t', comment='#', names=column_names, header=0)
    #Additionally separating "attributes" column on ";"
    attributes_split = df['attributes'].str.split(';', expand=True)
    new_column_names = ['featureID', 'FeatureName', 'LocusTag', 'Product', 'Dbxref', 'Gene']
    attributes_split.columns = new_column_names[:attributes_split.shape[1]]
    tag_values = {'featureID': 'ID=',
    'FeatureName': 'Name=',
    'LocusTag': 'locus_tag=',
    'Product': 'product=',
    'Dbxref': 'Dbxref=',
    'Gene': 'gene='}
    def assign_if_contains_specific_value(cell, specific_value):
        if cell and specific_value in cell:
            return cell
        return " "
    new_columns = {}
    for col in attributes_split.columns:
        specific_value = tag_values.get(col, '')
        attributes_split[col] = attributes_split[col].apply(assign_if_contains_specific_value, args=(specific_value,))
    
    df = pd.concat([df.drop(columns=['attributes']), attributes_split], axis=1)

    def __init__(self):
        pass

    def featcounting(self):
        #Analisis of features in GFF3 file (should work w/o FASTA file)
        cdscount = 0
        ncRNAcount = 0
        tRNAcount = 0
        oriCcount = 0
        regregcount = 0
        gapcount = 0
        for feat_type in self.df["type"]:
            if feat_type == "CDS":
                cdscount = cdscount + 1
            if feat_type == "ncRNA":
                ncRNAcount = ncRNAcount + 1
            if feat_type == "tRNA":
                tRNAcount = tRNAcount + 1
            if feat_type == "oriC":
                oriCcount = oriCcount + 1
            if feat_type == "regulatory_region":
                regregcount = regregcount + 1
            if feat_type == "gap":
                gapcount = gapcount + 1
        anno_gen_count = cdscount + ncRNAcount + tRNAcount + oriCcount + regregcount + gapcount
        global features
        features = {"CDS":cdscount, "ncRNA":ncRNAcount, "tRNA":tRNAcount, "oriC":oriCcount, "Regulatory Regions":regregcount, "Gaps":gapcount,
                    "Annotated Genes": anno_gen_count}
        return features

    def coganalyse(self):
        #Analisis of COGs in GFF3 file (should work w/o FASTA file)
        cog_count = 0
        COGtype_list = []
        COG_type_sum = []
        for attribute in self.df["Dbxref"]:
            if "Dbxref=COG" in attribute:
                cog_count = cog_count + 1
        for letter in "ABCDEFGHIJKLMNOPQRSTUWXYZ":
            for COG_type in self.df["Dbxref"]:
                if str("COG:" + str(letter) + ",") in COG_type:
                    COGtype_list.append("COG:" + str(letter))
        COGtype_count_d = dict(Counter(COGtype_list)) 
        cog_share = round((cog_count / features["CDS"]) * 100, 2)
        cog_data = {"COG count": cog_count, "COG share": cog_share}
        cog_data.update(COGtype_count_d)
        for COG_type in cog_data:
            if "COG:" in COG_type:
                COG_type_sum.append (cog_data[COG_type])
        for COG_type in cog_data:
            if "COG:" in COG_type:
                COG_type_share = (cog_data[COG_type] / sum(COG_type_sum)) * 100
                cog_data.update ({COG_type : "Count: " + str(cog_data[COG_type]) + " Share: " + str(round(COG_type_share, 2)) + " %"})
        return cog_data
    
    def creator_gff3(self):
        #Function to collect info from GFF3 file to create new MultiFASTA file 
        indexlist = (self.df.index[self.df['type']=='CDS'].tolist())
        with_cds_d = {}
        for element in indexlist:
            id_w_cds = self.df.iloc[element]['featureID']
            start_w_cds = self.df.iloc[element]['start']
            end_w_cds = self.df.iloc[element]['end']
            strang_w_cds = self.df.iloc[element]['strand']
            #print ("seqid:", seqid_w_cds, "start:", start_w_cds, "end:", end_w_cds)
            with_cds_d.update ({str(id_w_cds): str(start_w_cds) + " " + str(strang_w_cds) + " " +str(end_w_cds)})
        return with_cds_d
    
    def csv_creator(self):
        #CSV file creator
        self.df.to_csv("samplefile.csv", sep=" ", float_format=str, )
        sucmes = "seccess"

        return sucmes

































#    import csv
#        featureID = []
#        FeatureTyp =[]
#        Start = []
#        Stop =[]
#        Genname = []
#        Product = []
#        Strang = []
#        for attributes in self.df["attributes"]:
#            attrlist = str(attributes).split(";")
#            #print (attrlist)
#            for attrs in attrlist:
#                if "ID=" in str(attrs):
#                    featureID.append(attrs[3:])
#                if "Name=" in str(attrs):
#                    Genname.append(attrs[5:])
#                if "product" in str(attrs):
#                    Product.append(attrs[8:])
#               
#        for types in self.df["type"]:
#            FeatureTyp.append(str(types))
#        for starts in self.df["start"]:
#            Start.append(str(starts))
#        for stops in self.df["end"]:
#            Stop.append(str(stops))
#        for strangs in self.df["strand"]:
#            Strang.append(str(strangs))
#
#        with open('samplefile.csv', 'w', newline='') as csv_file:
#            csv_writer = csv.writer(csv_file)
#            csv_field = ["FeatureID", "FeatureTyp", "Start", "Stop", "Genname", "Product", "Strang"]
#
#            csv_writer.writerow(csv_field)
#            for n in range (0,10000):
#                try:
#                    csv_writer.writerow([featureID[n], FeatureTyp[n], Start[n], Stop[n], Genname[n], Product[n], Strang[n]])
#                except:
#                    break
#            sucmes = "success"
#            return sucmes






            

        








