from .urls import REQUEST_ADDRESS
from .urls import (TCGA_ESTIMATED_SIZE_FIELD,
                   TCGA_SUBMISSION_TIME_FIELD,
                   TCGA_TICKET_ID_FIELD,
                   TCGA_STATUS_CHECK_URL_FIELD)

import requests
import logging
import time
import os
import hashlib
import json

from appdirs import user_data_dir

PYTCGA_BASE_DIRECTORY = user_data_dir("pytcga", version="0.1")

def cache_data_dir():
    if not os.path.exists(PYTCGA_BASE_DIRECTORY):
        os.makedirs(PYTCGA_BASE_DIRECTORY)
    return PYTCGA_BASE_DIRECTORY

def check_if_exists_cached_file(output_file_name):
    return os.path.exists(
        os.path.join(cache_data_dir(), output_file_name)
    )

def tcga_request(disease,
                  center=None,
                  level=None,
                  platform=None,
                  platformType=None,
                  sample_list=None,
                  consolidateFiles='true',
                  flattenDir='true',
                  cache=True,
                  wait_time=30):

    # All disease codes must be upper-case.
    disease = disease.upper()

    filter_parameters = {
        'disease': disease,
        'center': center,
        'level': level,
        'platform': platform,
        'platformType': platformType,
        'sampleList': sample_list,
        'flattenDir': flattenDir,
        'consolidateFiles': consolidateFiles
    }

    # Save hash of request parameters
    request_id = hashlib.md5(
                        json.dumps(filter_parameters,
                                   sort_keys=True).encode('utf-8')
                    ).hexdigest()

    # Create an output tar file with that ID
    output_file_name = request_id + '.tar'

    # If using the cache, check if the file already exists
    if cache:
        archive_path = os.path.join(cache_data_dir(), output_file_name)

        if os.path.exists(archive_path):
            return archive_path

    (ticket_id, status_url) = create_tcga_request(filter_parameters)
    return check_and_retrieve_archive(
                                status_url,
                                archive_file_name=output_file_name,
                                wait_time=wait_time)


def create_tcga_filter_request(disease,
                              center=None,
                              level=None,
                              platform=None,
                              platformType=None,
                              sample_list=None,
                              consolidateFiles='true',
                              flattenDir='true'):
    """Create a web service request for TCGA

    Parameters
    ----------
    disease : str
        TCGA disease code {'LUAD', 'BRCA', 'BLAC',...}
    center : str
        TCGA center code {'BI', ...}
    level : str, optional
        Level 1, 2 or 3
    platform : str
        Platform {'Automated Mutation Calling', 'Curated Mutation Calling', ...}
    platformType : str
        Platform type {'Somatic Mutations', ...}
    sample_list : str, optional
        List of specific TCGA IDs to request
    consolidateFiles : str, optional
        Combine files of the same type, 'true' or 'false'
    flattenDir : str, optional
        Flatten the directory structure, 'true' or 'false'

    Raises
    ------
    ValueError
        A platform is required for the web service

    Returns
    -------
    (ticket_id, status_url) : (str, str)
        A pair of the ticket_id created and status url
    """
    if platform is None:
        raise ValueError(
            """Platform must be specfied."""
        )

    # All disease codes must be upper-case.
    disease = disease.upper()

    filter_parameters = {
        'disease': disease,
        'center': center,
        'level': level,
        'platform': platform,
        'platformType': platformType,
        'sampleList': sample_list,
        'flattenDir': flattenDir,
        'consolidateFiles': consolidateFiles
    }

    return create_tcga_request(filter_parameters)

class RequestError(Exception):
    def __init__(self, code, msg):
        self.code = code
        self.msg = msg

def create_tcga_request(filter_parameters):
    """Builds a webservice request from a dictionary of parameters

    Parameters
    ----------
    filter_parameters : dict
        Dictionary of filter parameters for TCGA web service request

    Returns
    -------
    (ticket_id, status_url) : (str, str)
        A pair of the ticket_id created and status url
    """
    response = requests.get(REQUEST_ADDRESS, params=filter_parameters)

    logging.debug("Request has status code {}".format(response.status_code))

    if response.status_code == 200:
        # Even if the servers sends a 200 code, it doesn't mean that
        # our request went through. For example, if you request data
        # from the wrong center, the server returns 200, but the page
        # will say: HTTP STATUS 204 - No Content...
        # The following checks for these non-std errors.
        rtext = response.text
        if rtext.startswith("<h2>HTTP STATUS"):
            error_code = rtext.split(' ')[2]
            error_msg = rtext.split(' - ')[1]
            raise RequestError(error_code, error_msg)

        # If all looks fine, go ahead and parse the JSON
        parsed_response = response.json()

        ticket_id = parsed_response[TCGA_TICKET_ID_FIELD]
        submission_time = parsed_response[TCGA_SUBMISSION_TIME_FIELD]
        estimated_size = parsed_response[TCGA_ESTIMATED_SIZE_FIELD]
        status_url = parsed_response[TCGA_STATUS_CHECK_URL_FIELD]

        logging.info(
            'Request received at {} for ticket {}, estimated size {}'.format(
                submission_time,
                ticket_id,
                estimated_size
            )
        )

        logging.info(
            'Track ticket {} at {}'.format(
                ticket_id,
                status_url
            )
        )
    else:
        raise ValueError('Request {} failed, \n{}'.format(response.url, response.text))

    return (ticket_id, status_url)

def retrieve_ticket_status(status_url):
    """Extract the current status at the tracking/status url

    Parameters
    ----------
    status_url : str
        Tracking URL for a submitted data request

    Returns
    -------
    job_status : str
        Current status {'OK', 'Accepted', ...}
    """
    tracking_response = requests.get(status_url)
    tracking_response_parsed = tracking_response.json()
    job_status = tracking_response_parsed['job-status']

    return job_status

def retrieve_archive(archive_url,
                    output_file_name,
                    block_size=1024):
    """Download the archive from given URL into the output file

    Parameters
    ----------
    archive_url : str
        TCGA URL to the created archive
    output_file_name : str
        Filename to save the archive

    Returns
    -------
    archive_path : str
        Full path to the downloaded archive
    """
    archive_path = os.path.join(cache_data_dir(), output_file_name)
    logging.info('Saving request to {}'.format(archive_path))

    with open(archive_path, 'wb') as archive:
        archive_response = requests.get(archive_url, stream=True)

        for block in archive_response.iter_content(block_size):
            archive.write(block)

    return archive_path

def check_and_retrieve_archive(status_url,
                               archive_file_name,
                               wait_time=None):
    """Checks the status URL from TCGA and attempts to download the archive

    Returns None if no archive exists and did not re-poll

    Parameters
    ----------
    status_url : str
        URL to check the status of the data request
    archive_file_name : str
        File name for archive from TCGA data request
    wait_time : int, optional
        Wait-time (in seconds) before polling service for data

    Returns
    -------
    archive_path : str
        Return a path to the downloaded archive
        Returns None if no archive exists and did not re-poll
    """
    job_status = retrieve_ticket_status(status_url)
    if job_status['status-message'] == 'OK':
        archive_url = job_status['archive-url']
        return retrieve_archive(archive_url, archive_file_name)

    if wait_time:
        while(True):
            time.sleep(wait_time)
            job_status = retrieve_ticket_status(status_url)

            if job_status['status-message'] == 'OK':
                archive_url = job_status['archive-url']
                archive_path = retrieve_archive(archive_url, archive_file_name)
                return archive_path

    return None
