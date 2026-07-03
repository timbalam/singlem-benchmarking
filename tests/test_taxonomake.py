
import os
import extern
import shutil
import pytest

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

@pytest.fixture
def end_to_end():
    def cleanup():
        try:
            os.remove(f"{path_to_data}/tmp/small_1.fq.gz")
        except FileNotFoundError:
            pass
        try:
            os.remove(f"{path_to_data}/tmp/small_2.fq.gz")
        except FileNotFoundError:
            pass
    
    cleanup()
    cmd = f"taxonomake {path_to_data}/community.toml"
    extern.run(cmd)
    yield
    cleanup()

def test_taxonomake(end_to_end):
    assert os.path.isfile(f"{path_to_data}/tmp/small_1.fq.gz")
    assert os.path.isfile(f"{path_to_data}/tmp/small_2.fq.gz")