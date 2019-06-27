import logging

from .vcf_tools import get_pop_AF

LOG = logging.getLogger(__name__)

UNINFORMATIVE       = 0
UPD_MATERNAL_ORIGIN = 1
UPD_PATERNAL_ORIGIN = 2
ANTI_UPD            = 3
PB_HOMOZYGOUS       = 4
PB_HETEROZYGOUS     = 5


def upd_site_call(gt_pb, gt_mo, gt_fa):
    """Call UPD informative sites
    
    input genotype coding: 0=HOM_REF, 1=HET, 3=HOM_ALT, 2=other

    Args:
        gt_pb (int): Genotype of proband as integer
        gt_mo (int): Genotype of mother as integer
        gt_fa (int): Genotype of father as integer
    
    Returns:
        site_info (int): One of the globals
    
    """

    opposite = {0:3, 3:0}
    
    # If any of the individuals are 'other' the site is UNINFORMATIVE
    if gt_pb == 2 or gt_mo == 2 or gt_fa == 2:
        return 0

    if gt_pb == 0 or gt_pb == 3:
        opp = opposite[gt_pb]
        if gt_mo == opp and gt_fa == opp:
            return UNINFORMATIVE
        if gt_mo == opp:
            return UPD_PATERNAL_ORIGIN
        if gt_fa == opp:
            return UPD_MATERNAL_ORIGIN
        return PB_HOMOZYGOUS
        
    if gt_pb == 1:
        if gt_mo != 1 and gt_fa != 1 and gt_mo == opposite[gt_fa]:
            return ANTI_UPD
        return PB_HETEROZYGOUS

    return UNINFORMATIVE


def get_UPD_informative_sites(vcf, csq_fields, proband, mother, father, min_af=0.05, 
                              af_tag='MAX_AF', min_gq=30):
    """Get UPD calls for each informative SNP above given pop freq
    
    Args:
        vcf (cyvcf2.VCF)
        min_af (float): Minimum allele frequency to consider SNP
        af_tag (str): Key to AF in annotation
        min_gq (int): Minimum GQ to consider variant
        csq_fields (list): describes VEP annotation
        sids (cyvcf2.samples): describes what position inds have in VCF
        proband (str): ID of proband in VCF
        mother (str): ID of mother in VCF
        father (str): ID of father in VCF
        
    Yields:
        site_calls (dict): A generator with dictionaries that describes the variant.

    """
    site_calls = []
    # Make UPD call for site and append to call list
    sids = vcf.samples
    proband_idx = sids.index(proband)
    mother_idx  = sids.index(mother)
    father_idx  = sids.index(father)
    
    nr_informative = 0
    
    for i,var in enumerate(vcf,1):
    
        # Raise error if multi-allelic site
        if len(var.ALT) > 1:
            raise SystemExit('ERROR: Split your variants!')

        # Skip non-SNPs
        if not var.is_snp:
            continue
        
        # Skip variants with population frequency < threshold
        if min_af > get_pop_AF(var, csq_fields, af_tag):
            continue
        
        # Skip variants where any individual has GQ < threshold
        if not all(gq >= min_gq for gq in var.gt_quals):
            continue

        gt = var.gt_types
        pos_call = upd_site_call(gt[proband_idx], gt[mother_idx], gt[father_idx])
        
        nr_informative += 1
        yield {'chrom':var.CHROM, 'pos':var.POS, 'call':pos_call}
    
    LOG.info("%s variants in vcf", i)
    LOG.info("%s informative variants found", nr_informative)

def call_regions(sites):
    """Yields called regions
    
    A putative_call is a collection of information from the variants in a upd region. 
    The call is started if there are variants that indicates UPD maternal or paternal origin. When 
    a anti upd site is encountered or the chromosome ends the call returned.
    
    Args:
        sites (iterable(dict))

    Yields:
        putative_call (dict): UPD informative region 
    """
    opposite = {UPD_MATERNAL_ORIGIN:UPD_PATERNAL_ORIGIN, UPD_PATERNAL_ORIGIN:UPD_MATERNAL_ORIGIN}
    calls = []
    putative_call = None
    prev = None
    last_seen_anti = {'chrom':'0', 'pos':0}
    
    for c in sites:
        if prev and (c['chrom'] != prev['chrom']):
            LOG.info("Chromosome %s checked", prev['chrom'])
        if not putative_call:
            # Create new putative call if UPD informative position
            if c["call"] in [UPD_MATERNAL_ORIGIN, UPD_PATERNAL_ORIGIN]:
                putative_call = {
                    'call':c['call'],
                    'chrom':c['chrom'],
                    'start_lo':last_seen_anti['pos'],
                    'start_hi':c['pos'],
                    'end_lo':c['pos'],
                    'end_hi':c['pos'],
                    'run_len': 1,
                    'opposites': 0,
                    'hom_sites': 1,
                    'het_sites': 0,
                    'tot': 1
                }

            # Save last positive which is definately not in UPD region
            if c["call"] == ANTI_UPD or not prev or c['chrom'] != prev['chrom']:
                last_seen_anti = {'chrom':c['chrom'], 'pos':c['pos']}

            prev = c
            continue

        # If an anti-UPD site or end of chromosome => close the putative call
        if c["call"] == ANTI_UPD or c['chrom'] != putative_call['chrom']:
            if c["call"] == ANTI_UPD:
                putative_call['end_hi'] = c['pos']-1
            else:
                putative_call['end_hi'] = prev['pos']
            
            yield putative_call
            
            last_seen_anti = {'chrom': c['chrom'], 'pos': c['pos']}      
            putative_call = None
            prev = c
            continue

        # If site call is same as putative call => extend the putative region
        if c['call'] == putative_call['call']:
            putative_call['end_lo'] = c['pos']
            putative_call['run_len'] += 1

        if c['call'] == PB_HOMOZYGOUS:
            putative_call['hom_sites'] += 1
        
        if c['call'] == PB_HETEROZYGOUS:
            putative_call['het_sites'] += 1
                    
        # If site call is opposite (maternal<->paternal), count it. (This pretty much never happens?)
        if c['call'] == opposite[putative_call['call']]:
            putative_call['opposites'] += 1

        # Count total number of SNPs in the call region
        putative_call['tot'] += 1

        prev = c

    LOG.info("Chromosome %s checked", prev['chrom'])

    if putative_call:
        yield putative_call