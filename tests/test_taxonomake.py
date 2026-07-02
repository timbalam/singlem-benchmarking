
import os
import extern

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

def test_taxonomake():
    cmd = "taxonomake {path_to_data}/community.toml --output {path_to_data}/tmp/test_taxonomake"
    extern.run(cmd)
