import urllib2
import json
import StringIO


def build_oncotator_search_string(
        chromosome_number,
        chromosome_start_pos,
        chromosome_end_pos,
        ref_allele,
        var_allele
        ):
    search_string = '{0}_{1}_{2}_{3}_{4}'.format(
        chromosome_number,
        chromosome_start_pos,
        chromosome_end_pos,
        ref_allele,
        var_allele
    )
    return search_string


def retrieve_oncotator_mutation_data_as_json(
        chromosome_number=None,
        chromosome_start_pos=None,
        chromosome_end_pos=None,
        ref_allele=None,
        var_allele=None,
        search_string=None
        ):
    """
    Parameters
    ----------
    chromosome_number: int
    chromosome_start_pos: int
    chromosome_end_pos: int
    ref_allele: str
    var_allele: str

    Returns
    -------
    oncotator_data: dict
    """
    if search_string is None:
        search_string = build_oncotator_search_string(
            chromosome_number,
            chromosome_start_pos,
            chromosome_end_pos,
            ref_allele,
            var_allele
        )
    page = retrieve_oncotator_mutation_data(search_string)
    return json.load(StringIO.StringIO(page))


def retrieve_oncotator_mutation_data(search_string_query, maxreadlength=100000000):
    base_url = 'http://www.broadinstitute.org/oncotator/mutation/{0}/'
    url_request_string = base_url.format(search_string_query)
    response = urllib2.urlopen(url_request_string, timeout=None)
    page = response.read(maxreadlength)
    return page
