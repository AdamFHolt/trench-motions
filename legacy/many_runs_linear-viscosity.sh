#!/bin/bash

ref_frame="nnr"
python compute_rates_misfit.py $ref_frame  7 1 23.5e6 1e-13 1 0
python compute_rates_misfit.py $ref_frame  7 1 58.7e6 1e-13 1 0
python compute_rates_misfit.py $ref_frame  7 1 93.9e6 1e-13 1 0

ref_frame="hs3"
python compute_rates_misfit.py $ref_frame  7 1 23.5e6 1e-13 1 0
python compute_rates_misfit.py $ref_frame  7 1 58.7e6 1e-13 1 0
python compute_rates_misfit.py $ref_frame  7 1 93.9e6 1e-13 1 0

ref_frame="sa"
python compute_rates_misfit.py $ref_frame  7 1 23.5e6 1e-13 1 0
python compute_rates_misfit.py $ref_frame  7 1 58.7e6 1e-13 1 0
python compute_rates_misfit.py $ref_frame  7 1 93.9e6 1e-13 1 0


# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 1 0
# python compute_rates_misfit.py  $ref_frame 7 1 23.5e6 1e-13 1 0
# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 1 0
# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 1 0
# python compute_rates_misfit.py  $ref_frame 4 1 23.5e6 1e-14 1 0
#python compute_rates_misfit.py  $ref_frame 4 1 23.5e6 1e-13 1 0
#python compute_rates_misfit.py  $ref_frame 4 1 23.5e6 1e-15 1 0

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 0.25

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 0.25

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 0.25

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 0.1

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 0.1

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 0.1

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 0.5

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 0.5

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 0.5

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 1
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 1
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 1

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 1
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 1
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 1

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 1
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 1
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 1


# ref_frame="sa"

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 0.25

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 0.25

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 0.25

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 0.1

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 0.1

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 0.1

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 0.5

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 0.5

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 0.5

# ref_frame="hs3"

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 0.25

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 0.25

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 0.25
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 0.25

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 0.1

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 0.1

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 0.1
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 0.1

# python compute_rates_misfit.py  $ref_frame 1 1 23.5e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 2 1 23.5e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 3 1 23.5e6 1e-13 0.5

# python compute_rates_misfit.py  $ref_frame 1 1 58.7e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 2 1 58.7e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 3 1 58.7e6 1e-13 0.5

# python compute_rates_misfit.py  $ref_frame 1 1 93.9e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 2 1 93.9e6 1e-13 0.5
# python compute_rates_misfit.py  $ref_frame 3 1 93.9e6 1e-13 0.5
