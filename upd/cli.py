import logging

import coloredlogs
import click

from pprint import pprint as pp

from upd.__version__ import __version__
from upd.vcf_tools import (parse_CSQ_header, get_vcf)
from upd.utils import (get_UPD_informative_sites, call_regions)
from upd.bed_utils import (output_filtered_regions)

LOG = logging.getLogger(__name__)

LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(__version__)
    ctx.exit()

@click.group()
@click.option('--vcf',
    type=click.Path(exists=True),
    required=True,
)
@click.option('--proband',
    help="ID of proband in VCF",
    required=True,
)
@click.option('--mother',
    help="ID of mother in VCF",
    required=True,
)
@click.option('--father',
    help="ID of father in VCF",
    required=True,
)
@click.option('--af-tag',
    help="Which field to use for population frequency filtering",
    default='MAX_AF',
    show_default=True
)
@click.option('--vep',
    help="If af-tag is in VEP annotation",
    is_flag=True,
)
@click.option('--min-af',
    help="Minimum SNP frequency",
    default=0.05,
    show_default=True
)
@click.option('--min-gq',
    help="Minimum GQ score",
    default=30,
    show_default=True
)
@click.option('--loglevel', 
    default='INFO', 
    type=click.Choice(LOG_LEVELS),
    help="Set the level of log output.", 
    show_default=True
)
@click.option('--version', 
    is_flag=True, 
    callback=print_version,
    expose_value=False, 
    is_eager=True
)
@click.pass_context
def cli(context, vcf, proband, mother, father, af_tag, vep, min_af, min_gq, loglevel):
    """Simple software to call UPD regions from germline exome/wgs trios"""
    coloredlogs.install(level=loglevel)
    LOG.info("Running upd version %s", __version__)

    context.obj = {}
    # Check if the given samples IDs exist in the VCF header
    try:
        vcf_reader = get_vcf(vcf, proband, mother, father)
    except Exception as err:
        LOG.warning(err)
        context.abort()

    csq_fields = None
    if vep:
        try:
            csq_fields = parse_CSQ_header(vcf_reader)
        except Exception as err:
            LOG.warning(err)
            context.abort()
        
        # Make sure the given VEP field exists
        if af_tag not in csq_fields:
            LOG.warning("The field %s does not exist in the VEP annotations", af_tag)
            LOG.info("Existing CSQ fields {}".format('|'.join(csq_fields)))
            context.abort()
    else:
        if not vcf_reader.contains(af_tag):
            LOG.warning("The field %s does not exist in the VCF", af_tag)
            context.abort()
    
    # Get all UPD informative sites into a list 
    context.obj['site_calls'] = get_UPD_informative_sites(
        vcf=vcf_reader,
        csq_fields=csq_fields,
        proband=proband,
        mother=mother,
        father=father,
        min_af=min_af,
        af_tag=af_tag,
        min_gq=min_gq
    )

@cli.command()
@click.option('--min-sites',
    help="Minimum UPD informative sites required to call a region",
    default=3,
    show_default=True
)
@click.option('--min-size',
    help="Minimum size (bp) required to call a region",
    default=1000,
    show_default=True
)
@click.option('-o','--out',
    help="Output bed file of all informative sites",
    type=click.Path(exists=False),
    default='-',
)
@click.pass_context
def regions(context, min_sites, min_size, out):
    """Call UPD regions"""
    # Make region calls
    calls = call_regions(context.obj['site_calls'])

    out_lines = output_filtered_regions(calls, min_sites, min_size)

    with click.open_file(out, 'w') as f:
        for line in out_lines:
            f.write(line+'\n')


@cli.command()
@click.option('-o','--out',
    help="Output bed file of all informative sites",
    type=click.Path(exists=False),
    default='-',
)
@click.pass_context
def sites(context, out):
    """Prints the sites that are informative for UPD"""
    site_type_name = [
        "UNINFORMATIVE", "UPD_MATERNAL_ORIGIN", "UPD_PATERNAL_ORIGIN", "ANTI_UPD", 
        "PB_HOMOZYGOUS", "PB_HETEROZYGOUS"
    ]
    
    with click.open_file(out, 'w') as f:
        for scall in context.obj['site_calls']:
            f.write("{}\t{}\t{}\t{}\n".format(
                scall['chrom'],
                scall['pos']-1,
                scall['pos'],
                site_type_name[scall['call']]
            ))
