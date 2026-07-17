
import os
import extern
import shutil
import pytest

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
snakefile = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "src/taxonomake/modules/Snakefile")

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
    yield
    cleanup()

def test_taxonomake(end_to_end):
    cmd = f"taxonomake {path_to_data}/community.yaml"
    extern.run(cmd)
    assert os.path.isfile(f"{path_to_data}/tmp/small_1.fq.gz")
    assert os.path.isfile(f"{path_to_data}/tmp/small_2.fq.gz")


@pytest.fixture
def end_to_end_gtdbtk():
    def cleanup():
        try:
            os.remove(f"{path_to_data}/tmp/small_1.fq.gz")
        except FileNotFoundError:
            pass
        try:
            os.remove(f"{path_to_data}/tmp/small_2.fq.gz")
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(f"{path_to_data}/tmp/genomes.gtdbtk")
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(f"{path_to_data}/tmp/gtdbtk_data")
        except FileNotFoundError:
            pass
    
    cleanup()
    yield
    cleanup()

def test_taxonomake_gtdbtk(end_to_end_gtdbtk):
    cmd = f"taxonomake {path_to_data}/community_gtdbtk.toml"
    extern.run(cmd)
    assert os.path.isfile(f"{path_to_data}/tmp/small_1.fq.gz")
    assert os.path.isfile(f"{path_to_data}/tmp/small_2.fq.gz")