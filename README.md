# kukulkan-
Convert Kaiju and Kraken2 outputs into BIOM format with taxonomy tables for downstream analysis.  Tool to transform Kaiju/Kraken2 classifications into valid BIOM files and taxonomic tables for Phyloseq &amp; vegan

###########################################################################################
KUKULKAN-tax v1.0 — #######################################################################
###########################################################################################
KUKULKAN-tax converts Kaiju classification reports into analysis-ready outputs for R, including a Phyloseq-compatible BIOM table and vegan-style abundance tables at multiple taxonomic ranks. It also exports a full taxonomy table mapped from NCBI’s rankedlineage.dmp.
__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

--What it does (at a glance)--

Parses NCBI rankedlineage.dmp to build a complete taxonomy look-up (Realm → Kingdom → Phylum → Class → Order → Family → Genus → Species → Name).
Reads multiple Kaiju report files (one per sample), counts only classified reads (status C) per TaxID, and combines them into an OTU/feature table.
Builds vegan-style tables for each taxonomic level (rows = samples; columns = taxa at that level).
Creates a BIOM v1.0 JSON file with observation metadata (taxonomy lineage) that Phyloseq can import directly.
__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

--Inputs

1. Kaiju reports list (-i/--input)
A plain-text file where each line is an absolute path to a Kaiju per-sample report (TSV).
  -Expected 3 columns per report: Status, Read_name, TaxID (tab-separated).
  -Only rows where Status == "C" (classified) are counted.

Sample name rule: the script uses the file stem and takes everything before the first underscore _ as the sample ID.
Example: file S01_kaiju_output.txt → sample name S01.

2. NCBI taxonomy file (-t/--taxonomy-file)
Absolute path to rankedlineage.dmp (from NCBI new_taxdump).
You may provide a virus-only subset (e.g., rankedlineage_viruses.dmp) if you’ve pre-filtered it.

Output folder (-o/--output)
Absolute path to an existing directory where results will be written.
Important: The script concatenates paths as strings, so end your path with a trailing / (e.g., /path/to/out/).

Flag for unclassified handling (--na-uniq-id)
If present, taxa that are NA at a given level get a unique placeholder per TaxID (Unc.<TaxID>) in vegan tables, letting you keep different “unclassified” taxa separated.
Note: In the current script, this flag is marked as required by the argument parser; you must include it to run (presence = enabled).

__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
Outputs (written to --output)

tax_table.tsv — Full taxonomy table (index = TaxID; columns = Realm … Species, Name).
otu_table.tsv — Counts matrix (rows = TaxID, columns = samples) with integer counts of classified reads.
vegan_<LEVEL>.tsv — One table per taxonomic level (Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species, Name):
Rows = samples, columns = taxa names at that level.
If --na-uniq-id is used, missing ranks are labeled Unc.<TaxID> to keep them distinct.
otu_biom.tsv — BIOM v1.0 JSON file (type "OTU table") embedding lineage metadata for each TaxID.

You can validate it with:

bash
biom validate-table -i otu_biom.tsv
__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

Installation and dependencies

  Python 3.8+ recommended

  Packages: pandas, biom-format

Quick setup (conda + pip)

bash
# create a clean environment (optional but recommended)
conda create -y -n kukulkan python=3.11
conda activate kukulkan

# install dependencies
pip install pandas biom-format
__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

Getting the NCBI taxonomy file

Download the current rankedlineage.dmp from the NCBI new_taxdump:

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/rankedlineage.dmp


You may pre-filter it to keep only viruses if desired (producing, e.g., rankedlineage_viruses.dmp).

Expected format (pipe-separated):
TaxID | Name | Species | Genus | Family | Order | Class | Phylum | Kingdom | Superkingdom

KUKULKAN-tax internally reverses these columns to a top-down order: Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species, Name.

__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

Preparing the Kaiju paths file

Create a text file (e.g., kaiju_paths.txt) with one absolute path per line:


/data/kaiju/S01_kaiju_output.txt
/data/kaiju/S02_kaiju_output.txt
/data/kaiju/S03_kaiju_output.txt

Make sure the filenames follow the SAMPLE_*.txt pattern so the script can extract SAMPLE correctly.
__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

Usage
python kukulkan_tax.v1.py \
  --input /abs/path/kaiju_paths.txt \
  --taxonomy-file /abs/path/rankedlineage.dmp \
  --output /abs/path/outdir/ \
  --na-uniq-id


Notes:

The output directory must already exist and must end with /.

Because the current parser marks --na-uniq-id as required, you must include it; its presence enables unique Unc.<TaxID> labels for missing ranks in vegan tables. (If you prefer the old behavior—omitting the flag to keep a single shared NA/Unclassified—set the argument to not be required in the script.)

__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

How it works (details)

Load taxonomy
The script reads rankedlineage.dmp, trims whitespace, replaces blanks with NA, sets TaxID as index, and relabels columns to ['Realm','Kingdom','Phylum','Class','Order','Family','Genus','Species','Name']. It writes this as tax_table.tsv and keeps an in-memory dictionary for fast lookups. 

Build OTU table from Kaiju
For each report:

Keep only rows with Status == 'C'.

Count reads by TaxID (per sample).

Outer-merge counts across samples → otu_table.tsv (missing = 0). 

Generate vegan tables
For each level (Realm…Name), the script:

Maps each non-zero TaxID in a sample to its name at that level.

If --na-uniq-id is present, replaces NA with Unc.<TaxID> so “unclassified” bins remain distinct.

Sums counts per resulting taxon name (columns) to produce a sample-by-taxon table for that level. 

Create BIOM (v1.0 JSON)

Intersects the OTU TaxIDs with available taxonomy rows to keep them aligned.

Adds observation_metadata = {'taxonomy': [Realm…Name]} for each TaxID.

Exports a BIOM v1.0 JSON (saved as otu_biom.tsv). 


__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
Practical tips & conventions

Sample naming: The sample ID is everything before the first _ in each Kaiju filename. Keep this consistent.

Unclassified handling:

With --na-uniq-id: “unclassified at level” taxa are kept separate per TaxID (e.g., Unc.123456 vs Unc.7891011).

Without it (requires editing the script to make the flag optional): all such entries collapse into a single NA/Unclassified column per level.

Phyloseq import: Use biomformat::read_biom() (R) then coerce to a phyloseq object, or import via phyloseq::import_biom().

Validation: Always run biom validate-table -i otu_biom.tsv if you plan to share downstream.
__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


Example end to end  
# 1) environment
conda create -y -n kukulkan python=3.11
conda activate kukulkan
pip install pandas biom-format

# 2) taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/rankedlineage.dmp

# 3) list of Kaiju outputs (absolute paths)
printf "/data/kaiju/S01_kaiju_output.txt\n/data/kaiju/S02_kaiju_output.txt\n" > /data/kaiju/kaiju_paths.txt

# 4) run (note the trailing slash in --output)
python kukulkan_tax.v1.py \
  --input /data/kaiju/kaiju_paths.txt \
  --taxonomy-file /data/rankedlineage.dmp \
  --output /data/kukulkan_out/ \
  --na-uniq-id

# 5) validate BIOM
biom validate-table -i /data/kukulkan_out/otu_biom.tsv

__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
Troubleshooting

“No such file or directory” — Ensure all paths in --input are absolute and exist; ensure --output exists and ends with /.

“argument --na-uniq-id is required” — The current script requires this flag. Include --na-uniq-id or edit the parser to make it optional.

Empty vegan tables — Check that Kaiju reports contain Status == 'C' rows and that TaxIDs exist in your rankedlineage.dmp subset.

License & Attribution

KUKULKAN-tax v1.0 script and this README description are based on the provided source file. Please cite appropriately when using in publications.
_________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
Why use BIOM and taxonomy tables?

Most downstream microbiome/virome analyses in R and Python rely on standardized data structures rather than raw Kaiju reports. KUKULKAN-tax bridges that gap by producing two complementary output formats:

1. BIOM object (otu_biom.tsv)

  The BIOM format (Biological Observation Matrix) is a community standard supported by many tools.

  In R, BIOM files can be imported directly into Phyloseq, which allows:

    Alpha diversity metrics (Shannon, Simpson, Chao1, etc.)

    Beta diversity metrics (Bray–Curtis, Jaccard, weighted/unweighted UniFrac, etc.)

    Ordinations (PCoA, NMDS), PERMANOVA tests, rarefaction curves

    Taxonomic barplots, heatmaps, and advanced ecological visualizations

  This makes the BIOM output ideal if you want a single integrated object (counts + taxonomy) for interactive analysis and plotting.


2. Vegan-style abundance tables (vegan_<LEVEL>.tsv)

  These are flat TSV tables with samples in rows and taxa in columns, provided separately for each taxonomic level (Realm → Species).

  The format is simple and flexible, which makes it compatible with many R and Python packages, such as:

    vegan (community ecology: ordination, PERMANOVA, rarefaction)

    MaAsLin2 / MaAsLin3 (multivariable association testing)

    SparCC / SpiecEasi (microbial/viral network inference)

    Custom pipelines for correlation, regression, or machine learning

  These tables are especially useful if you want fine-grained control over which taxonomic level to analyze, or if you plan to test associations with clinical/phenotypic metadata.


3. Why this matters

Raw Kaiju outputs are not analysis-ready: they contain per-read classifications without aggregation or taxonomic hierarchy. For many users (especially those new to microbial/viral bioinformatics), manually converting those outputs into valid OTU tables, taxonomy tables, and Phyloseq objects is error-prone and time-consuming.


KUKULKAN-tax automates this entire step:

  Ensures taxonomic lineage consistency via NCBI rankedlineage.dmp.

  Produces validated BIOM files and level-specific abundance tables.

  Gives you immediately usable inputs for the most common R packages used in microbiome/virome research.

In short:

Use the BIOM file if you want to work in Phyloseq and leverage its ecosystem for diversity metrics and visualization.

Use the vegan-style tables if you want to run association or network analyses (MaAsLin2, SparCC, SpiecEasi, vegan).




