# workflow
# 

import nipype

# Build the confound design matrix
# use feat.py

# use fsf file to find the in_file
# 

# Regress the confounds out of the timeseries
#confregress = pe.MapNode(fsl.FilterRegressor(filter_all=True),
#                         iterfield=["in_file", "design_file", "mask"],
#                         name="confregress")
#