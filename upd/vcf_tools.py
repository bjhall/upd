
def get_header_desc(reader, header_id):
    """Get description field of an header field ID
    
    Args:
        reader (cyvcf2.VCF)
        header_id (str)
    
    Returns:
        str: Information from a vcf header description
    """
    for rec in reader.header_iter():
        d = rec.info()
        if d.get('ID') == header_id:
            return d.get('Description')

    return None

def parse_CSQ_header(reader):
    """Parse the order of VEP fields from the vcf header
    
    Args:
        reader (cyvcf2.VCF)

    Returns:
        csq_format (list(str)): A list with the VEP header
    """
    csq_str = get_header_desc(reader, "CSQ")
    
    if not csq_str:
        raise ValueError("CSQ header field missing. The VCF need to be annotated with VEP")
        
    _, csq_format_str = csq_str.split('Format: ')
    csq_format_str = csq_format_str.rstrip('"')
    csq_format = csq_format_str.split('|')
    return csq_format

def get_pop_AF(vep_data, vep_fields, AF_str):
    """Extract population frequency from VEP annotations.
    
    Args:
        vep_data (str): VEP annotation from vcf INFO field
        vep_fields (list): Description of VEP annotation
        AF_str (str): Name of AF field to parse
    
    Returns:
        freq (float): The annotated frequency, returns 0 if no data
    """

    first_vep_str = vep_data.split(',')[0]
    data = first_vep_str.split('|')

    for i in range(len(data)):
        if vep_fields[i] == AF_str:
            return float(data[i] or 0)

    return 0
