# ---------------------------------------------------------------------------------------------------------- #
#                                                  Librerías                                                 #
# ---------------------------------------------------------------------------------------------------------- #

import pandas as pd
from collections import Counter
import sys, argparse
from biom.table import Table
from datetime import datetime
import json
from pathlib import Path

# ---------------------------------------------------------------------------------------------------------- #
#                                                 ARGUMENTOS                                                 #
# ---------------------------------------------------------------------------------------------------------- #

parser = argparse.ArgumentParser(prog = "KUKULCAN tax (version 1)", description= "This program converts kaiju report data to phyloseq data.")

parser.add_argument('-i', '--input', required=True, help="File with paths of kaiju reports.")
parser.add_argument('-t', '--taxonomy-file', required=True, help="rankedlineage.dmp file")
parser.add_argument('-o', '--output', required=True, help="Output folder for results.")
parser.add_argument('--na-uniq-id', action='store_true', required=True, help="Create unique ID for each NA present (default: False).")

args = parser.parse_args()
infile = args.input
taxfile = args.taxonomy_file
outdir = args.output
unique_id = args.na_uniq_id

# ---------------------------------------------------------------------------------------------------------- #
#                                                  PRINCIPAL                                                 #
# ---------------------------------------------------------------------------------------------------------- #

# Create taxa table from ncbi data

tax_ncbi = pd.read_table(taxfile, index_col=0, header=None, sep='|').iloc[:, :-1]
tax_ncbi = tax_ncbi.map(lambda x: x.strip().replace(' ', '_') if isinstance(x, str) else x)
tax_ncbi.replace(r'^\s*$', 'NA', regex=True, inplace=True)

tax_ncbi.index.name = 'TaxID' 
tax_ncbi.columns = ['Name', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Realm']

tax_ncbi = tax_ncbi[tax_ncbi.columns[::-1]]
outfile_taxa = outdir + 'tax_table.tsv'
tax_ncbi.to_csv(outfile_taxa, sep='\t', index=True)

tax_dict = tax_ncbi.to_dict(orient='index')

# Create otu table with counts for each organism classified in each sample

with open(infile, "r") as file:
    kaiju_paths = [line.strip() for line in file if line.strip()]

for i, filepath in enumerate(kaiju_paths):
    filename = Path(filepath).stem
    sample = filename.split("_")[0]
    sample_data = pd.read_table(filepath, header=None, names=['Status', 'Read_name', 'TaxID'], sep="\t")
    sample_classified = sample_data[sample_data['Status'] == 'C']
    sample_counts = Counter(sample_classified['TaxID'])
    sample_table = pd.DataFrame(list(sample_counts.items()), columns=['TaxID', f"{sample}"])
    if i == 0:
        otu_table = sample_table 
    else:
        otu_table = pd.merge(otu_table, sample_table, on='TaxID', how='outer')

otu_table = otu_table.fillna(0).astype(int)
otu_table = otu_table.set_index('TaxID')
outfile_otu = outdir + 'otu_table.tsv'
otu_table.to_csv(outfile_otu, sep='\t', index=True)

# Create vegan files

tax_levels = list(tax_ncbi.columns)

for levels in tax_levels:
    first = True
    for col, data in otu_table.items():
        otu_keys = data[data != 0].index.tolist()

        if unique_id:
            sample_taxonomy = {
                taxid: {
                    taxa: (name if name != 'NA' else 'Unc.' + str(taxid)) 
                    for taxa, name in tax_dict[taxid].items() if taxa == levels
                }
                for taxid in otu_keys if taxid in tax_dict
            }
        else:
            sample_taxonomy = {
                taxid: {taxa: name for taxa, name in tax_dict[taxid].items() if taxa == levels}
                for taxid in otu_keys if taxid in tax_dict
            }

        sample_tax_data = pd.DataFrame.from_dict(sample_taxonomy, orient='index')
        sample_tax_counts = Counter(sample_tax_data[levels])
        sample_tax_table = pd.DataFrame([sample_tax_counts], index=[str(data.name)])

        if first:
            tax_table = sample_tax_table
            first = False
        else:
            tax_table = pd.concat([tax_table, sample_tax_table], axis=0, join='outer')

    tax_table = tax_table.fillna(0).astype(int)
    tax_table.index.name = 'Sample'
    outfile_vegan = outdir + 'vegan_' + levels + ".tsv"
    tax_table.to_csv(outfile_vegan, sep='\t', index=True)

# Asegurar taxIDs comunes
common_taxids = otu_table.index.intersection(tax_ncbi.index)
otu = otu_table.loc[common_taxids]
tax = tax_ncbi.loc[common_taxids]

# Construir metadata compatible con BIOM v1.0
metadata = []
for _, row in tax.iterrows():
    lineage = [str(i) if pd.notna(i) else "NA" for i in row]
    metadata.append({'taxonomy': lineage})

# Crear tabla biom
biom_table = Table(
    otu.values,
    observation_ids=list(otu.index),
    sample_ids=list(otu.columns),
    observation_metadata=metadata,
    type="OTU table"
)

# Generar JSON de la tabla
biom_json = biom_table.to_json("kaiju2phyloseq")

# Convertir a diccionario para asegurar estructura
biom_dict = json.loads(biom_json)

# Agregar campos obligatorios del formato 1.0
biom_dict["format"] = "Biological Observation Matrix 1.0.0"
biom_dict["format_url"] = "http://biom-format.org"
biom_dict["generated_by"] = "kaiju2phyloseq pipeline"
biom_dict["date"] = datetime.now().isoformat()
biom_dict["id"] = "otu_table"

outfile_biom = outdir + 'otu_biom.tsv'

# Guardar como archivo BIOM válido
with open(outfile_biom, 'w') as f:
    json.dump(biom_dict, f, indent=4)

print("✅ Archivo .biom compatible con BIOM v1.0 generado exitosamente.")
