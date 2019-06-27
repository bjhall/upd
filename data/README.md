test.vcf.gz is a randomly generated example vcf to be used for testing. Execute the following:

```bash
upd --vcf data/test.vcf.gz --proband TEST_PROBAND --mother TEST_MOTHER --father TEST_FATHER --vep regions
```

Expected output:

```
15      22958118        101910531       ORIGIN=PATERNAL;TYPE=HETERODISOMY;LOW_SIZE=78952412;INF_SITES=138;SNPS=2032;HET_HOM=889/1006;OPP_SITES=0;START_LOW=20170126;END_HIGH=102516586;HIGH_SIZE=82346460
```