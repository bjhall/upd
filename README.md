## UPD.py

Simple software to call UPD regions from germline exome/wgs trios.

## Documentation

### Installation

The only dependency is cyvcf2. Only tested on python 2.7...


### Input VCF requirements

* Trio (proband-mother-father) VCF is required (i.e. all three individuals in the same VCF).
* GQ must be present in genotype fields.
* Must be annotated with VEP and have some sort of population frequency field in the CSQ annotation (default: MAX_AF)


### Running

```bash
upd.py input.vcf.gz --proband PB_ID --mother MOTHER_ID --father FATHER_ID --out upd_regions.bed
```

Where PB_ID/MOTHER_ID/FATHER_ID are the sample IDs from the vcf header.

#### Optional parameters
--min_af (DEFAULT: 0.05)
Minimum population frequency required to include the SNP in the analysis.

--vep_af (DEFAULT: MAX_AF)
Specifies which CSQ (typically from VEP) field to be used when selecting SNPs. Getting population frequencies from other INFO fields are currenly not supported.

--min_gq (DEFAULT: 30)
Specifies the minimum GQ required to include a variant in the analysis. All three individuals' must have a GQ larged than or equal to this.

--min_sites (DEFAULT: 3)
Minimum number of consecutive UPD sites needed to call an UPD region.

--min_size (DEFAULT: 1000)
Minimum number of base pairs between first and last UPD site in a region required to call it.

--out_sites FILENAME
If given a BED file with all informative sites used to call regions will be written.


### Output files

#### Called regions (--out)
BED file of all called regions fulfilling the --max_sites and --max_size filter criteria, including some annotation.

15      22958103        101910550       ORIGIN=PATERNAL;TYPE=HETERODISOMY;LOW_SIZE=78952446;INF_SITES=182;SNPS=2896;HET_HOM=1275/1440;OPP_SITES=0;START_LOW=20170150;END_HIGH=102516492;HIGH_SIZE=82346342


Field | Description
-------------------
ORIGIN | The parental origin of the UPD region
TYPE | Indicates whether region is heterodisomic or isodisomic/deletion. See caveats
LOW_SIZE | Lower bound of the estimated size of the called region. Distance between first and last UPD site of the region
INF_SITES | Number of UPD informative sites in the called regions
SNPS | Total number of SNPs in the region
HET_HOM | Number of heterozygous and homozygous sites in the region. Used to indicate if region has LOH
OPP_SITES | Number of sites in the region that have the opposite parental origin. So far I've never seen any other number than 0...
START_LOW | Earliest possible start positon of the UPD region (position of last anti-UPD site prior to the region)
END_HIGH | Last possible end position of the UPD region (position of first anti-UPD site prior to the region)
HIGH_SIZE | Upper bound of the estimated size of the called region. Distance between the surrounding anti-UPD sites.

#### Caveats
The isodisomic/heterodisomic calling depends on estimating if there is a run of homozygousity in the called region. This is called using presence of heterozygous sites within the call, and should only be seen as a rough estimate for smaller calls. For more statistically sound detection of this, use e.g. bcftools roh to detect regions of homozygozity and combine the results.

Isodisomies cannot be differentiated from deletions using this software. Combine the calls with SV/CNV calls to determine if there is an overlapping deletion. If overlapping deletion, the "UPD" call can be used to detemine which parental allele was deleted.




