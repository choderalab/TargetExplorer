import urllib2
import json
import StringIO

def retrieve_oncotator_mutation_data_as_json(search_string_query):
    page = retrieve_oncotator_mutation_data(search_string_query)
    return json.load(StringIO.StringIO(page))

def retrieve_oncotator_mutation_data(search_string_query, maxreadlength=100000000):
    base_url = 'http://www.broadinstitute.org/oncotator/mutation/{0}/'
    url_request_string = base_url.format(search_string_query)
    # url_request_string = base_url + urllib.quote_plus(search_string_query) + '&force=yes&format=xml'
    response = urllib2.urlopen(url_request_string)
    page = response.read(maxreadlength)
    return page