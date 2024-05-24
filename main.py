#!/usr/bin/env python3
import argparse
import datei as dt
import datei2 as dt2





def main():
    #ArgParser
    parser = argparse.ArgumentParser(description="Save the file path and run the reader script.")
    parser.add_argument("-f", nargs=1, type=str, help='The path to the FASTA', default=None)
    parser.add_argument("-g", nargs=1,type=str, help='The path to the GFF3', default=None)
    parser.add_argument("-b", nargs=2,type=str, help='The path to the both FASTA and GFFF3', default=None,)
    args = parser.parse_args()
    #print (args)

    #FNA FIle path management
    path_fna = args.f
    if args.f == None and args.b != None:
        path_fna = args.b[0]
    elif args.b == None and args.f !=None:
        path_fna = args.f
    else:
        print ("No FASTA file provided")
    path_fna = str(path_fna).strip()
    path_fna = path_fna.strip("[']")

    #GFF3 File path management
    path_gff3 = args.g
    if args.g == None and args.b != None:
        path_gff3 = args.b[1]
    elif args.b == None and args.g != None:
        path_gff3 = args.g
    else:
        print ("No GFF3 file provided")
    path_gff3 = str(path_gff3).strip()
    path_gff3 = path_gff3.strip("[']")

    #Path saver
    with open('config.txt', 'w') as file:
        file.write(f'path_fna={path_fna}\n')
        file.write(f'path_gff3={path_gff3}\n')
    
    #main stuff
        gff3file = dt2.GFF3seq()
        featcount = gff3file.featcounting()
        print ("Features:")
        for f in featcount:
            print (f, ":", featcount[f])
        cog_data = gff3file.coganalyse ()
        print ("COGS:", cog_data["COG count"], "(", cog_data["COG share"], "%)")
        print ("COG Categories:")
        for COG_type in cog_data:
            if "COG:" in COG_type:
                print (COG_type, ":", cog_data[COG_type])
        testiter = gff3file.creator_gff3()
        print (len(testiter))
        csv_maker = gff3file.csv_creator()
        print (csv_maker)
        

        fnafile = dt.FASTAseq()
        seq_w_feat = fnafile.creator_fasta()
        
        start_cod_dict, all_cod_dict = fnafile.codone_usage()
        print ("Codone Usage Analysis:")
        print ("Start Codone Usage Analysis:")
        for k,v in start_cod_dict.items():
            print ("Start codone: ", k, "Count:", v)
        print ("All-Codone Usage Analysis:")
        for k,v in all_cod_dict.items():
            if len(k) == 3:
                print ("Codone: ", k, "Count: ", v)
        plot = fnafile.visualiser_codones()
        #print (seq_w_feat)
        #test gitlab change
        

            

if __name__ == "__main__":
    main()
     


