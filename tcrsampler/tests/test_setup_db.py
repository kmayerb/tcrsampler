import pytest
from tcrsampler import setup_db

def test_install_all_next_gen(dry_run = True):
	""" with dry_run set to true, test setup_db routine runs without error"""
	setup_db.install_all_next_gen(dry_run = dry_run)

def test_install_nextgen_data_to_db_KeyError_download_from(download_file = "new_nextgen_chains_mouse_A.tsv", 
									download_from = "wrong_from", 
									dry_run = True):
	""" with dry_run set to true, test setup_db routine runs without error"""
	with pytest.raises(KeyError):
		setup_db.install_nextgen_data_to_db(download_file = download_file, 
											download_from = download_from , 
											dry_run = dry_run)

def test_install_nextgen_data_to_db_KeyError_filename(download_file = "new_nextgen_chains_hamster_A.tsv", 
									download_from = "dropbox", 
									dry_run = True):
	""" with dry_run set to true, test setup_db routine runs without error"""
	with pytest.raises(KeyError):
		setup_db.install_nextgen_data_to_db(download_file = download_file, 
											download_from = download_from , 
											dry_run = dry_run)

def test_install_nextgen_data_to_db_TypeError(download_file = 1, 
									download_from = "dropbox", 
									dry_run = True):
	""" with dry_run set to true, test setup_db routine runs without error"""
	with pytest.raises(TypeError) as e:
		setup_db.install_nextgen_data_to_db(download_file = download_file, 
										download_from = download_from , 
										dry_run = dry_run)
		assert e.type is TypeError
		assert e.value.args[0].startswith('The <download_from> arg must be a string')


def test_install_nextgen_dry_run_true(download_file = 'ruggiero_mouse_sampler.zip', 
									  download_from = "dropbox", 
									  dry_run = False):
	""" with dry_run set to true, test setup_db routine runs without error"""
	setup_db.install_nextgen_data_to_db(download_file = download_file, 
										download_from = download_from , 
										dry_run = dry_run)

