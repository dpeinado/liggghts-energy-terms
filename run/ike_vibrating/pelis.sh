#!/bin/bash

python ~/liggghts-energy-terms/pyPost/pyVTK.py files/dump.Hooke_en_0.85_cof_0.5 vtks/hooke_en_0.85_cof_0.5_
python ~/liggghts-energy-terms/pyPost/pyVTK.py files/dump.Hooke_en_0.30_cof_1.0 vtks/hooke_en_0.30_cof_1.0_
python ~/liggghts-energy-terms/pyPost/pyVTK.py files/dump.Hertz_en_0.85_cof_0.5 vtks/hertz_en_0.85_cof_0.5_
python ~/liggghts-energy-terms/pyPost/pyVTK.py files/dump.Hertz_en_0.30_cof_1.0 vtks/hertz_en_0.30_cof_1.0_
