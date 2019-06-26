## UPD.py

Simple software to call UPD regions from germline exome/wgs trios.

## Documentation

### Installation

`pip install --editable .` or `python setup.py install`

The only dependencies are cyvcf2, click and coloredlogs. If cyvcf2 does not install please try to install it 
manually with [bioconda](https://anaconda.org/bioconda/cyvcf2) and then run the command above.


### Input VCF requirements

* Trio (proband-mother-father) VCF is required (i.e. all three individuals in the same VCF).
* GQ must be present in genotype fields.
* Must be annotated with some sort of population frequency. This can also be in the CSQ annotation (default: MAX_AF), if so use `--vep`


### Running

```bash
upd --vcf input.vcf.gz --proband PB_ID --mother MOTHER_ID --father FATHER_ID regions --out upd_regions.bed
```

Where PB_ID/MOTHER_ID/FATHER_ID are the sample IDs from the vcf header.

#### Optional parameters
Command |Parameter | Description
------- |--------- | -----------
base | **--min-af (DEFAULT: 0.05)** | Minimum population frequency required to include the SNP in the analysis.
base | **--af-tag (DEFAULT: MAX_AF)** | Specifies which frequency field to be used when selecting SNPs.
base | **--min-gq (DEFAULT: 30)** | Specifies the minimum GQ required to include a variant in the analysis. All three individuals' must have a GQ larged than or equal to this.
base | **--vep (flag)** | If given, search the CSQ field for `af-tag`
base | **--vep (flag)** | If given, search the CSQ field for `af-tag`
regions | **--min-sites (DEFAULT: 3)** | Minimum number of consecutive UPD sites needed to call an UPD region.
regions | **--min-size (DEFAULT: 1000)** | Minimum number of base pairs between first and last UPD site in a region required to call it.
regions/sites | **--out (DEFAULT: stdout)** | If the results should be printed to a file


### Output

#### Called regions (udp regions)
BED file of all called regions fulfilling the --min-sites and --min-size filter criteria, including some annotation.

Example call:
```
15      22958103        101910550       ORIGIN=PATERNAL;TYPE=HETERODISOMY;LOW_SIZE=78952446;INF_SITES=182;SNPS=2896;HET_HOM=1275/1440;OPP_SITES=0;START_LOW=20170150;END_HIGH=102516492;HIGH_SIZE=82346342
```

Field | Description
----- | -----------
ORIGIN | The parental origin of the UPD region
TYPE | Indicates whether region is heterodisomic or isodisomic/deletion. See caveats
LOW_SIZE | Lower bound of the estimated size of the called region. Distance between first and last UPD site of the region
INF_SITES | Number of UPD informative sites in the called regions
SNPS | Total number of SNPs in the region
HET_HOM | Number of heterozygous and homozygous sites in the region. Used to indicate if region has LOH
OPP_SITES | Number of sites in the region that have the opposite parental origin. So no other number than 0 have been observed...
START_LOW | Earliest possible start positon of the UPD region (position of last anti-UPD site prior to the region)
END_HIGH | Last possible end position of the UPD region (position of first anti-UPD site prior to the region)
HIGH_SIZE | Upper bound of the estimated size of the called region. Distance between the surrounding anti-UPD sites.

#### Informative sites (upd sites)
Fourth field indicates what type of field:

Type | Description
---- | -----------
ANTI_UPD | Site where proband has inherited alleles from both parents (e.g. AA x BB = AB)
UPD_MATERNAL_ORIGIN | Site where both alleles are inherited from the mother, or a fathers allele was deleted. (e.g. AA x BB = AA)
UPD_PATERNAL_ORIGIN | Site where both alleles are inherited from the father, or a mothers allele was deleted. (e.g. AA x BB = BB)
PB_HETEROZYGOUS | Heterozygous site in the proband (used only to call hetero/isodisomy)
PB_HOMOZYGOUS | Homozygous site in the proband (used only to call hetero/isodisomy)
UNINFORMATIVE | Various sites excluded from the analysis. Ignore these...


### Caveats
The isodisomic/heterodisomic calling depends on estimating if there is a run of homozygousity in the called region. This is called using presence of heterozygous sites within the call, and should only be seen as a rough estimate for smaller calls. For more statistically sound detection of this, use e.g. bcftools roh to detect regions of homozygozity and combine the results.

Isodisomies cannot be differentiated from deletions using this software. Combine the calls with SV/CNV calls to determine if there is an overlapping deletion. If overlapping deletion, the "UPD" call can be used to detemine which parental allele was deleted.




