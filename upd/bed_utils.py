
def output_filtered_regions(calls, min_sites=3, min_size=1000, iso_het_pct=0.01):
    """Takes called regions, filters them, and yields annotated BED lines

    Args:
        calls (iterable): An iterable with UPD regions
        min_sites (int): Minimum UPD informative sites required to call a region
        min_size (int): Minimum size (bp) required to call a region
        iso_het_pct(float): Quota between het-sites and hom-sties
    Yields:
        out_line (str): A formated string with the information about a region
    """
    out_line = ("{}\t{}\t{}\tORIGIN={};TYPE={};LOW_SIZE={};INF_SITES={};SNPS={};HET_HOM={}/{}"
               ";OPP_SITES={};START_LOW={};END_HIGH={};HIGH_SIZE={}")

    for rcall in calls:
        if not rcall['run_len'] > min_sites:
            continue
        call_str = {2: 'PATERNAL', 1: 'MATERNAL'}
        lo_size = rcall['end_lo'] - rcall['start_hi']
        hi_size = rcall['end_hi'] - rcall['start_lo']

        # Rough estimation if heterodisomy or homodisomy/deletion
        het_pct = float(rcall['het_sites']) / float((rcall['het_sites']+rcall['hom_sites']))
        upd_type = "HETERODISOMY"
        if het_pct < iso_het_pct:
            upd_type = "HOMODISOMY/DELETION"

        if lo_size < min_size:
            continue

        yield out_line.format(
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
            )
