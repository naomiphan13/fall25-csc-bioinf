# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import os

def data_dir():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    week3_root = os.path.abspath(os.path.join(current_dir, "../../.."))
    data_dir = os.path.join(week3_root, "data")
    return data_dir
