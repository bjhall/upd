#!/usr/bin/env python2
from cyvcf2 import VCF
import argparse

# Constants
UNINFORMATIVE       = 0
UPD_MATERNAL_ORIGIN = 1
UPD_PATERNAL_ORIGIN = 2
ANTI_UPD            = 3
PB_HOMOZYGOUS       = 4
PB_HETEROZYGOUS     = 5
site_type_name = ["UNINFORMATIVE", "UPD_MATERNAL_ORIGIN", "UPD_PATERNAL_ORIGIN", "ANTI_UPD", "PB_HOMOZYGOUS", "PB_HETEROZYGOUS"]

def main():

    opt = parse_arguments()

    vcf_reader = VCF(opt.vcf)

    # Check if the given samples IDs exist in the VCF header
    sids = vcf_reader.samples
    if not all(elem in sids for elem in [opt.proband, opt.mother, opt.father]):
        raise SystemExit("ERROR: At least one of the given sample IDs do not exist in the VCF header")

    # Parse the VEP CSQ header
    csq_fields = parse_CSQ_header(vcf_reader)

    # Make sure the given VEP field exists
    if opt.vep_af not in csq_fields:
        raise SystemExit("ERROR: The field {} does not exist in the VEP annotations".format(opt.vep_af))

    # Get all UPD informative sites into a list 
    site_calls = get_UPD_informative_sites(vcf_reader, opt, csq_fields, sids)

    # Make region calls
    calls = call_regions(site_calls)

    # Output BED file with results
    output_filtered_regions(opt.out, calls, opt.min_sites, opt.min_size)

    # Output individual informative sites if requested
    if opt.out_sites:
        output_informative_sites(opt.out_sites, site_calls)

def get_UPD_informative_sites(vcf, opt, csq_fields, sids):
    """Get UPD calls for each informative SNP above given pop freq"""
    site_calls = []
    for var in vcf:
    
        # Raise error if multi-allelic site
        if len(var.ALT) > 1:
            raise SystemExit('ERROR: Split your variants!')

        # Skip non-SNPs
        if not var.is_snp:
            continue
        
        # Skip variants with population frequency < threshold
        if opt.min_af > get_pop_AF(var.INFO["CSQ"], csq_fields, opt.vep_af):
            continue
        
        # Skip variants where any individual has GQ < threshold
        if not all(gq >= opt.min_gq for gq in var.gt_quals):
            continue

        # Make UPD call for site and append to call list
        proband_idx = sids.index(opt.proband)
        mother_idx  = sids.index(opt.mother)
        father_idx  = sids.index(opt.father)

        gt = var.gt_types
        pos_call = upd_site_call(gt[proband_idx], gt[mother_idx], gt[father_idx])
        site_calls.append( {'chrom':var.CHROM, 'pos':var.POS, 'call':pos_call} )
    
    return site_calls


def call_regions(sites):
    """Takes a list of UPD informative sites and generates called regions"""
    opposite = {UPD_MATERNAL_ORIGIN:UPD_PATERNAL_ORIGIN, UPD_PATERNAL_ORIGIN:UPD_MATERNAL_ORIGIN}
    calls = []
    putative_call, prev = None, None
    last_seen_anti = 0

    for c in sites:
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

        else:

            # If an anti-UPD site or end of chromosome => close the putative call
            if c["call"] == ANTI_UPD or c['chrom'] != putative_call['chrom']:
                if c["call"] == ANTI_UPD:
                    putative_call['end_hi'] = c['pos']-1
                else:
                    putative_call['end_hi'] = prev['pos']
                calls.append(putative_call)
                last_seen_anti = {'chrom': c['chrom'], 'pos': c['pos']}      
                putative_call = None

            else:
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

    return calls


def output_filtered_regions(outfile, calls, min_sites, min_size):
    """Takes called regions, filters them, and outputs annotated BED file"""
    with open(outfile, 'w') as out:
    
        for rcall in calls:
            if rcall['run_len'] > min_sites:
                call_str = {2: 'PATERNAL', 1: 'MATERNAL'}
                lo_size = rcall['end_lo'] - rcall['start_hi']
                hi_size = rcall['end_hi'] - rcall['start_lo']

                # Rough estimation if heterodisomy or homodisomy/deletion
                het_pct = float(rcall['het_sites']) / float((rcall['het_sites']+rcall['hom_sites']))
                upd_type = "HETERODISOMY"
                if het_pct < 0.01:
                    upd_type = "HOMODISOMY/DELETION"
                if lo_size >= min_size:
                    out.write("{}\t{}\t{}\tORIGIN={};TYPE={};LOW_SIZE={};INF_SITES={};SNPS={};HET_HOM={}/{};OPP_SITES={};START_LOW={};END_HIGH={};HIGH_SIZE={}\n".format(
                        rcall['chrom'],
                        rcall['start_hi']-1,
                        rcall['end_lo'],
                        call_str[rcall['call']],
                        upd_type,
                        lo_size,
                        rcall['run_len'],
                        rcall['tot'],
                        rcall['het_sites'],
                        rcall['hom_sites'],
                        rcall['opposites'],
                        rcall['start_lo'],
                        rcall['end_hi'],
                        hi_size
                    ))



def output_informative_sites(outfile, calls):
    """Outputs all called informative sites to a BED file"""
    with open(outfile, 'w') as out:
    
        for scall in calls:
            
            out.write("{}\t{}\t{}\t{}\n".format(
                scall['chrom'],
                scall['pos']-1,
                scall['pos'],
                site_type_name[ scall['call'] ]
            ))
                    

def upd_site_call(gt_pb, gt_mo, gt_fa):
    """ Call UPD informative sites, input genotype coding: 0=HOM_REF, 1=HET, 3=HOM_ALT, 2=other """

    opposite = {0:3, 3:0}
    
    if gt_pb == 2 or gt_mo == 2 or gt_fa == 2:
        return UNINFORMATIVE

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


def get_header_desc(reader, header_id):
    """ Get description field of an header field ID """
    for rec in reader.header_iter():
        d = rec.info()
        if d.get('ID') == header_id:
            return d.get('Description')

    return None


def parse_CSQ_header(reader):
    """Parse the order of VEP fields from the vcf header"""
    csq_str = get_header_desc(reader, "CSQ")

    if csq_str:
        _, csq_format_str = csq_str.split('Format: ')
        csq_format_str = csq_format_str.rstrip('"')
        csq_format = csq_format_str.split('|')
        return csq_format

    else:
        raise ValueError("CSQ header field missing. The VCF need to be annotated with VEP")

    
def get_pop_AF(vep_data, vep_fields, AF_str):
    """Extract population frequency from VEP annotations, returns 0 if no data."""

    first_vep_str = vep_data.split(',')[0]
    data = first_vep_str.split('|')

    for i in range(len(data)):
        if vep_fields[i] == AF_str:
            return float(data[i] or 0)

    return 0


def parse_arguments():
    parser = argparse.ArgumentParser(description='Aggregate vcfs from different variant callers.')
    parser.add_argument('vcf', metavar='vcf', type=str, help='input vcf file')
    parser.add_argument('--proband', type=str, help='ID of proband in VCF', required=True)
    parser.add_argument('--mother', type=str, help='ID of mother in VCF', required=True)
    parser.add_argument('--father', type=str, help='ID of father in VCF', required=True)
    parser.add_argument('--vep_af', type=str, default='MAX_AF', help='Which VEP fields to use for poulation frequency filtering')
    parser.add_argument('--min_af', type=float, default=0.05, help='Minimum SNP frequency')
    parser.add_argument('--min_gq', type=int, default=30, help='Minimum GQ score')
    parser.add_argument('--min_sites', type=int, default=3, help='Minimum UPD informative sites required to call a region')
    parser.add_argument('--min_size', type=int, default=1000, help='Minimum size (bp) required to call a region')
    parser.add_argument('--out', type=str, required=True, help='Output bed file of called regions')
    parser.add_argument('--out_sites', type=str, default=None, help='Output bed file of all informative sites')

    opt = parser.parse_args()
    return opt


if __name__ == '__main__':
    main()
