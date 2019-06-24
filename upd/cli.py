import logging
# import coloredlogs
import click

from pprint import pprint as pp

from cyvcf2 import VCF

from upd.__version__ import __version__
from upd.vcf_tools import parse_CSQ_header
from upd.utils import (get_UPD_informative_sites, call_regions)
from upd.bed_utils import (output_filtered_regions, output_informative_sites)

LOG = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
LOG.addHandler(handler)

LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(__version__)
    ctx.exit()

@click.command()
@click.argument('vcf')
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
@click.option('--vep-af',
    help="Which VEP fields to use for population frequency filtering",
    default='MAX_AF',
    show_default=True
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
    help="Output bed file of called regions",
    type=click.Path(exists=False),
    required=True
)
@click.option('--out-sites',
    help="Output bed file of all informative sites",
    type=click.Path(exists=False),
)
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
@click.option('--loglevel', default='INFO', type=click.Choice(LOG_LEVELS),
              help="Set the level of log output.", show_default=True)
@click.pass_context
def cli(context, vcf, proband, mother, father, vep_af, min_af, min_gq, min_sites, min_size, out, 
        out_sites, loglevel):
    """Annotate str variants with str status"""
    # coloredlogs.install(level=loglevel)
    LOG.setLevel(loglevel)
    LOG.info("Running upd version %s", __version__)
    
    vcf_reader = VCF(vcf)
    
    # Check if the given samples IDs exist in the VCF header
    sids = vcf_reader.samples
    if not all(elem in sids for elem in [proband, mother, father]):
        context.abort("At least one of the given sample IDs do not exist in the VCF header")
    
    # Parse the VEP CSQ header
    try:
        csq_fields = parse_CSQ_header(vcf_reader)
    except ValueError as err:
        context.abort(err)
    
    # Make sure the given VEP field exists
    if vep_af not in csq_fields:
        LOG.warning("The field {} does not exist in the VEP annotations".format(vep_af))
        LOG.info("Existing CSQ fields {}".format('|'.join(csq_fields)))
        context.abort()
    
    # Get all UPD informative sites into a list 
    site_calls = get_UPD_informative_sites(
        vcf=vcf_reader, 
        csq_fields=csq_fields, 
        sids=sids, 
        proband=proband, 
        mother=mother, 
        father=father, 
        min_af=min_af, 
        vep_af=vep_af, 
        min_gq=min_gq
    )
        
    # Make region calls
    calls = call_regions(site_calls)

    # Output BED file with results
    output_filtered_regions(out, calls, min_sites, min_size)
    
    # Output individual informative sites if requested
    if out_sites:
        output_informative_sites(out_sites, site_calls)
    
    