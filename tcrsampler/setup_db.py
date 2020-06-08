import os
import sys

def install_nextgen_data_to_db(download_file, download_from = "dropbox", dry_run = False):
    """
    Function installs next-gen

    Parameters
    ----------
    download_file : string
        "new_nextgen_chains_mouse_A.tsv"
    download_from : string 
        'dropbox' is only currently available option
    dry_run : bool 
        if True, download commands are printed but no executed
    Returns
    -------
    curl_url_cmd : string
        the command for installing to tcrdist/db/alphabeta_db.tsv_files/*

    """
    if not isinstance(download_from, str):
        raise TypeError('The <download_from> arg must be a string')
    
    download_sources = ["dropbox"]
    if download_from not in download_sources:
        raise KeyError(f"the <download_from> arg must be be one of the following {download_sources}")

    if download_from is "dropbox":
        address = { "new_nextgen_chains_mouse_A.tsv" : 'https://www.dropbox.com/s/pkpr6p97eworn3q/new_nextgen_chains_mouse_A.tsv?dl=1',
                    "new_nextgen_chains_mouse_B.tsv" : 'https://www.dropbox.com/s/sxgvrj25mnzr20s/new_nextgen_chains_mouse_B.tsv?dl=1',
                    "new_nextgen_chains_human_A.tsv" : 'https://www.dropbox.com/s/9p43c7tscf46dat/new_nextgen_chains_human_A.tsv?dl=1',
                    "new_nextgen_chains_human_B.tsv" : 'https://www.dropbox.com/s/83xbk7jp4xd3pr0/new_nextgen_chains_human_B.tsv?dl=1',
                    "new_nextgen_chains_human_G.tsv" : 'https://www.dropbox.com/s/41w8yl38nr4ey32/new_nextgen_chains_human_A.tsv?dl=1',
                    "new_nextgen_chains_human_D.tsv" : 'https://www.dropbox.com/s/8ysciqrcywdsryp/new_nextgen_chains_human_B.tsv?dl=1'}

        if download_file not in address.keys():
            raise KeyError("download_file must be in {}".format(" ".join(map(str,address.keys))))

    # Where the file is to be installed
    path_file = os.path.realpath(__file__)
    path = os.path.dirname(path_file)
    install_path = os.path.join(path, "db", download_file) #<----- ALPHA BETA destination

    def generate_curl(filename, download_link):
        return('curl -o {} {} -L'.format(filename, download_link))

    curl_url_cmd = generate_curl(install_path, address[download_file])
    sys.stdout.write("RUNNING: {}\n".format(curl_url_cmd) )
    
    if dry_run is False:
        os.system(curl_url_cmd)
    return(curl_url_cmd)


def install_all_next_gen(dry_run = False):
    for fn in [ "new_nextgen_chains_mouse_A.tsv",
                "new_nextgen_chains_mouse_B.tsv" ,
                "new_nextgen_chains_human_A.tsv",
                "new_nextgen_chains_human_B.tsv",
                "new_nextgen_chains_human_G.tsv", 
                "new_nextgen_chains_human_D.tsv" ]:

        install_nextgen_data_to_db(download_file = fn, download_from = 'dropbox', dry_run = dry_run )

