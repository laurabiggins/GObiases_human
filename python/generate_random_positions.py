import random
from collections import defaultdict
#import pandas

# default output filename if one is not supplied
output_filename = "closest_genes.txt"
# number of genes to generate per chromosome
total_number_of_genes = 200
# file with lengths of chromosomes
chr_list = "D:/GObiases_human/data-raw/chr_list_noMT_HS.txt"
#gtf_file = "D:/GObiases_human/data-raw/Homo_sapiens.GRCh38.98_gene_info.txt"
#gtf_file = "D:/GObiases_human/data-raw/gene_info_head.txt"
gtf_file = "D:/Homo_sapiens.GRCh38.99.gtf"
#gtf_file = "D:/GObiases_human/data-raw/head.gtf"
biotype = "any"

total_number_of_genes = int(total_number_of_genes)
  
# dictionary of random positions, keys are chr, values are lists of positions
random_positions = {}
#genome_sizes = pandas.read_csv("D:/GObiases_human/data-raw/chr_list_noMT.csv", header=None, names = ("chr", "length"))

with open(chr_list) as f:
    genome_sizes = f.readlines()

# make this a float so that we can do float maths to determine how many genes to generate per chr
total_genome_size = float(0)

for line in genome_sizes:
    row = line.split("\t")
    # add the chr_length
    chr_length = int(row[1])
    total_genome_size += chr_length
    print(f'the length of chr {row[0]} is {chr_length:,}')


# total = f'Amount including tax is £{pounds + (tax_rate*0.01*pounds):,}'
print(f'Total genome size = {total_genome_size:,}')    


# we want the total number of genes to be approximately the number of genes specified, 
# it doesn't have to be exact (the rounding may mean we don't get the exact value)
actual_total_number_of_genes = int(0)
    
# generating the random positions
for line in genome_sizes:

    row = line.split("\t")
    
    chr = row[0]
    chr_length = int(row[1])
    
    # the number of random positions we'll generate for this chr
    number_of_genes = int(round((chr_length/total_genome_size)*total_number_of_genes))
   
    actual_total_number_of_genes += number_of_genes
    
    # create a list of random positions for each chr
    random_pos = []
    
    for x in range(number_of_genes):
        random_pos.append(random.randint(1,chr_length))
    
    random_positions[chr] = random_pos
    print ("%s random positions generated for chr %s" % (len(random_pos), chr))

string = "total number of random positions generated: %s "
print(string %(actual_total_number_of_genes))    

# load in gtf file
with open(gtf_file) as f:
    gtf_file = f.readlines()
    

# a dictionary of lists of lists to store gene info for each chromosome
# had to change from a list of tuples as these are immutable and I wanted to 
# include the distance from the between the random position and the gene
# keys are chr, values are gene info
gene_info = defaultdict(list)
       
for line in gtf_file:
    
    line = line.rstrip()
    if not line.startswith("#"):
    
        split_line=line.split("\t")
        
        # check whether the line contains a gene
        if (split_line[2] == "gene"):
            
            # parse the detail column that contains the name
            details = split_line[8].split(";")
            
            if("gene_biotype" in details[4]):               
                gene_biotype = details[4].replace("gene_biotype ", "")
                gene_biotype = gene_biotype.replace('"', "")
                gene_biotype = gene_biotype.replace(' ', "")
            else:
                print("gene biotype not found in gtf file")   
            
            # if biotype has been selected
            if((biotype=="protein_coding" and gene_biotype=="protein_coding") or biotype=="any"):
  
                if("gene_id" in details[0]):            
                    gene_id = details[0].replace("gene_id ", "")
                    gene_id = gene_id.replace('"', "")
                else:
                    print("gene id not found in gtf file")
                
                if("gene_name" in details[2]):   
                    gene_name = details[2].replace("gene_name ", "")
                    gene_name = gene_name.replace('"', "")
                    gene_name = gene_name.replace(' ', "")
                    #print(f'gene name = {gene_name}')
                else:
                    print("gene name not found in gtf file")
                   

                # load the gene info into a tuple and add to the list
                # in the gtf column 0 is chr, 3 is start, 4 is end
                chr = split_line[0]
                info = [gene_id, gene_name,chr,int(split_line[3]),int(split_line[4]),gene_biotype]
                
                #print(f'Info for chr{chr} is {info}') 
                
                gene_info[chr].append(info) 
                
            
# create a list of tuples to contain the closest genes
closest_genes = []


# go through each chr
for chr in random_positions:
        
    for random_pos in random_positions[chr]:
    
        closest_pos = 1000000000
        closest_gene = [] # this is a list so that it can contain the gene, the position of the gene and
        # the random position
        
        for gene in gene_info[chr]:

            # first check whether the random position is within a gene
            if (gene[3] < random_pos and gene[4] > random_pos):
                            
                closest_pos=0
                closest_gene = list(gene)
                
                # exit the for loop as we've found an overlapping gene
                break
                
            diff1 = abs(int(gene[3] - random_pos))
            diff2 = abs(int(gene[4] - random_pos))
            
            min_diff = min([diff1,diff2])
            
            # if smaller than the previous minimum, replace it  
            if(min_diff <= closest_pos):
                closest_pos = min_diff
                closest_gene = list(gene)
     
        closest_gene.append(closest_pos)
        closest_gene.append(random_pos)
        
        closest_genes.append(closest_gene)
          
header_line = "\t".join(["gene_id","gene_name","chromosome","start","end","biotype","distance", "random_pos"])

# write out the closest genes
with open(output_filename, "w") as f:
    
    f.write(header_line)
    f.write("\n")
    f.write("\n".join(["\t".join([str(g) for g in closest_gene]) for closest_gene in closest_genes]))

