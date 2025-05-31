#!/usr/bin/env python3

"""Assignment was to process three tiers of evidence to annotate FASTA protein sequences from files provided in class: 
   1. Hidden Markov Model results 
   2. BLAST results from a custom MySQL database
   3. TMHMM 2.0 results"""

import re # for pulling polypeptide id's from faa file and TMHMM results file
import pymysql.cursors # for processing BLAST results (TIER 2)


faa_products = {}

# create empty dictionary of all expected polypeptide ids from fasta file with empty list values to hold first matching product

for line in open("prodigal2fasta.nostars.faa"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        m = re.search(r"1_[0-9]+", line)
        faa_id = m.group()
        faa_products[faa_id] = []


# TIER 1: Process HMM Results (tab-delimited file) and add products meeting e-value < 1e-50 criteria

for line in open("hmmscan.htab"):
    linelist = line.split("\t")
    
    # remove _polypeptide ending from qry_id to match faa_id keys
    clean_id = linelist[5].replace("_polypeptide", "")
    product = linelist[15]
    evalue = linelist[19]
    # filter for results with evalues < 1e-50. Only add first product to the value list for each polypeptide id
    if float(evalue) < 10**-50 and len(faa_products[clean_id]) == 0:
        faa_products[clean_id].append(product)

       
# TIER 2: Process BLAST results from MySQL to find any products that meet the evalue <1e-50 criteria in BLAST but not in HMM

# Get MySQL username, password, and host from user
conn = pymysql.connect(user = input('Enter MySQL user: '), password = input('Enter user password: '), host = input('Enter host (eg. localhost): '), database='annot')
with conn:
    with conn.cursor() as cursor:

    # Filter MySQL BLAST results based on evalue < 1e-50 criteria
        qry = ("""
            SELECT qry_id, evalue, product
            FROM blast WHERE evalue < %s and product IS NOT NULL
            """)
        cursor.execute(qry, ("1e-50",))
        
        for (qry_id, evalue, product) in cursor:
            if len(faa_products[qry_id]) == 0:
                faa_products[qry_id].append(str(product))
        
        conn.commit()
        cursor.close()

# TIER 3: Process TMHMM results (tab-delimited file) to annotate any polypeptides that have transmembrane helices and aren't already annotated

for line in open("prodigal2fasta.nostars.tmhmm.short"):
    linelist = line.split("\t")
    clean_id = linelist[0].replace("_polypeptide", "")

    # extract number of helices from result file
    m = re.search(r"[0-9]+",linelist[2])
    helices = int(m.group())

    if len(faa_products[clean_id]) == 0 and helices > 0:
        faa_products[clean_id].append("Putative transmembrane protein")


# Write faa ids and annotated products to an output file
output_file = open("mpitt_annotated_final.txt", "w")

for faa_id, product in faa_products.items():
    if product == []:
        product.append("hypothetical protein")

    # make sure that there's only one annotated product per faa_id
    if len(product) == 1:
        output_file.write("{0}\t{1}\n".format(faa_id, "".join(product)))
    else:
        print("Improperly annotated polypeptide {0}".format(faa_id))

output_file.close()
