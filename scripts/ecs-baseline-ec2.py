# coding: utf-8
#!/usr/bin/env python
# coding: utf-8

# (c) British Crown Copyright, the Met Office.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following
# conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials
#       provided with the distribution.
#     * Neither the name of the Met Office nor the names of its
#       contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# RCS Id $Id: ecs-baseline-ec2.py,v 1.113 2020/06/04 13:24:34 hadmw Exp $"
#
# Modified from Sherwood et al., 2020:
# An Assessment of Earth's Climate Sensitivity Using Multiple Lines of Evidence
# doi: https://doi.org/10.1029/2019RG000678
# Code doi: http://doi.org/10.5281/zenodo.3945276

import numpy as np
import pickle
from joblib import dump, load
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import datetime
from emergent_constraints_lambda_parameters import *
from bellouin_ringberg_stats_fix_200311_ERFaer import *
import sys

print (sys.argv[0],"RCS ID=","$Id: ecs-baseline-ec2.py,v 1.113 2020/06/04 13:24:34 hadmw Exp $")

print ("args=", sys.argv)
calc_id = sys.argv[1]
outpath = sys.argv[2] + '/'
inpath  = sys.argv[3] + '/'
ref_paper = sys.argv[4].lower()

if 'sherwood' in ref_paper:
  # total cloud feedback, 0.45±0.33 (Sherwood et al., 2020)
  l_cld_mu = 0.45
  l_cld_sd = 0.33
else: # 'myers' in sys.argv[4].lower()
  # total cloud feedback, 0.27±0.25 (Myers et al., 2021)
  l_cld_mu = 0.27
  l_cld_sd = 0.25

out_dir = os.path.join(outpath, ref_paper, calc_id)
if not os.path.exists(out_dir):
  os.makedirs(out_dir)

# set numpy random seed to allow reproducible results
np.random.seed(42)

# in most cases calc_id is prior_id
prior_id = calc_id

# version with process EC deactivated for everything except UL_PROC
activate_process_ec = False

#set other defaults for cases where calc_id=prior_id
inflation=0

small_sample_pattern_effect=False
classic_hist_pdf=False
calc_tecs=False
fat_tails=False
lambda_ec=True
pattern_effect_dependence2=False
reduce_dt_dn_uncert=False
no_fhist_uncert=False
reduce_fhist_uncert=False
no_pattern_effect=False
half_pattern_effect=False
cmip5_pattern_effect=False

hist_forcing_version='Bellouin_2020_eqn_8'
#hist_forcing_version='AR5_extended'
#hist_forcing_version='Bellouin_2020_constrained'

plot_priors=0
plot_ecs_pdf=1

plot_prior=1
plot_posterior=1
plot_process_likelihood=0
plot_paleo_likelihood=0

no_process_bu=0
no_process_ec=1
no_hist=0
no_paleo_cold=0
no_paleo_hot=0

upper_medium_sample=0
medium_sample=0
prior_sample=0
tiny_sample=0
small_sample=0
quick_and_dirty=0
weight_prior=0
plot_frequency=1
pattern_effect_dependence=False
alternative_hist_dt=False
alternative_pi_period=False
period_1850_2005_2015=False
period_1750_2018=False

l_posterior_label=calc_id
s_posterior_label=calc_id

#exceptions to this are dealt with here

if calc_id == 'UL_INF':
  prior_id='UL'
  inflation=1

if calc_id == 'US1_INF':
  prior_id='US1'
  inflation=1

if calc_id == 'UL_HIST':
  prior_id='UL'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_HIST':
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_PRIOR':
  prior_id='US1'
  no_hist=1
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1

if calc_id == 'US1_HIST_NOPAT_NO_FHIST_UNCERT':
  no_fhist_uncert=True
  no_pattern_effect=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_HIST_NOPAT_REDUCE_FHIST_UNCERT':
  reduce_fhist_uncert=True
  no_pattern_effect=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_HIST_NOPAT':
  no_pattern_effect=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_HIST_SMALLSAMPLEPAT':
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  small_sample_pattern_effect=True

if calc_id == 'UL_SMALLSAMPLEPAT':
  prior_id='UL'
  small_sample_pattern_effect=True

if calc_id == 'US1_HIST_HALFPAT':
  half_pattern_effect=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_HIST_CMIP5PAT':
  cmip5_pattern_effect=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_HIST_ALT1':
  # Table 4.3 2nd Row
  no_pattern_effect=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  alternative_hist_dt=True
  alternative_DT_o=0.96  # update 19/3/20
  alternative_DT_sd=0.14/1.64

if calc_id == 'US1_HIST_ALT3':
  # Table 4.3 4th Row
  # Use ='AR5_extended' aerosol forcing
  hist_forcing_version='AR5_extended'
  no_pattern_effect=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_HIST_ALT4':
  # Table 4.3 4th Row
  # Use Bellouin constrained aerosol forcing
  hist_forcing_version='Bellouin_2020_constrained'
  no_pattern_effect=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_HIST_NOPAT_REDUCE_DTDN_UNCERT':
  no_pattern_effect=True
  reduce_dt_dn_uncert=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'US1_HIST_ALT2':
  # Table 4.3 3rd Row
  no_pattern_effect=True
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  alternative_hist_dt=True
  alternative_DT_o=1.02 # Updated 19/3/20
  alternative_DT_sd=0.12/1.64
  alternative_pi_period=True

if calc_id == 'UL_PALEO':
  prior_id='UL'
  no_hist=1
  no_paleo_cold=0
  no_paleo_hot=0
  no_process_bu=1
  no_process_ec=1
  no_hist=1
  plot_paleo_likelihood=1

if calc_id == 'US1_PALEO':
  prior_id='US1'
  no_hist=1
  no_paleo_cold=0
  no_paleo_hot=0
  no_process_bu=1
  no_process_ec=1
  no_hist=1
  plot_paleo_likelihood=1

if calc_id == 'US1_PALEOCOLD':
  prior_id='US1'
  no_hist=1
  no_paleo_cold=0
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  no_hist=1
  plot_paleo_likelihood=1

if calc_id == 'US1_PALEOHOT':
  prior_id='US1'
  no_hist=1
  no_paleo_cold=1
  no_paleo_hot=0
  no_process_bu=1
  no_process_ec=1
  no_hist=1
  plot_paleo_likelihood=1

if calc_id == 'UL_PROC':
  prior_id='UL'
  no_hist=1
  no_paleo_cold=1
  no_paleo_hot=1
  no_hist=1
  plot_process_likelihood=1

if calc_id == 'UL_PROCBU':
  prior_id='UL'
  no_hist=1
  no_paleo_cold=1
  no_paleo_hot=1
  no_hist=1
  no_process_bu=0
  no_process_ec=1
  plot_process_likelihood=1

if calc_id == 'US1_PROCBU':
  prior_id='US1'
  no_hist=1
  no_paleo_cold=1
  no_paleo_hot=1
  no_hist=1
  no_process_bu=0
  no_process_ec=1
  plot_process_likelihood=1

if calc_id == 'UL_PROCEC':
  prior_id='UL'
  no_hist=0
  no_paleo_cold=0
  no_paleo_hot=0
  no_hist=0
  no_process_bu=0
  no_process_ec=0
  plot_process_likelihood=1

if calc_id == 'UL_NOHIST':
  prior_id='UL'
  no_hist=1

if calc_id == 'UL_NOHIST_TECS':
  calc_tecs=True
  prior_id='UL'
  no_hist=1

if calc_id == 'UL_NOPALEO':
  prior_id='UL'
  no_paleo_cold=1
  no_paleo_hot=1

if calc_id == 'ULI_NOPALEO':
  no_paleo_cold=1
  no_paleo_hot=1
  prior_id = 'ULI'
  restrict_li_range=1
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  small_sample=0
  medium_sample=1
  tiny_sample=0
  ntop=50000
  weight_prior=0
  plot_frequency=100

if calc_id == 'UL_NOPALEOCOLD':
  prior_id='UL'
  no_paleo_cold=1

if calc_id == 'UL_NOPALEOCOLD_TECS':
  calc_tecs=True
  prior_id='UL'
  no_paleo_cold=1

if calc_id == 'UL_NOPALEOHOT':
  prior_id='UL'
  no_paleo_hot=1

if calc_id == 'UL_NOPALEOHOT_TECS':
  calc_tecs=True
  prior_id='UL'
  no_paleo_hot=1

if calc_id == 'US1_TECS':
  calc_tecs=True
  prior_id='US1'

if calc_id == 'US1_FAT':
  fat_tails=True
  prior_id='US1'

if calc_id == 'US1_NOPROCBU':
  prior_id='US1'
  no_process_bu=1

if calc_id == 'UL_NOPROCBU':
  prior_id='UL'
  no_process_bu=1

if calc_id == 'UL_NOPROCEC':
  prior_id='UL'
  no_process_ec=1

if prior_id == 'PROCBU':
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  small_sample=1
  tiny_sample=0
  ntop=10000
  weight_prior=0
  plot_frequency=ntop/10

if prior_id == 'UL':
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  small_sample=1
  tiny_sample=0
  ntop=10000
  weight_prior=0
  plot_frequency=ntop/10

if calc_id == 'UL_ALL_PROCEC':
  prior_id='UL'
  no_process_ec=0
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  small_sample=1
  tiny_sample=0
  ntop=10000
  weight_prior=0
  plot_frequency=ntop/10

if calc_id == 'UL_TECS':
  calc_tecs=True
  prior_id='UL'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  small_sample=1
  tiny_sample=0
  ntop=10000
  weight_prior=0
  plot_frequency=ntop/10

if calc_id == 'UL_FAT':
  fat_tails=True
  prior_id='UL'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  small_sample=1
  tiny_sample=0
  ntop=10000
  weight_prior=0
  plot_frequency=ntop/10

if calc_id == 'PLOT_ERF':
  prior_id='UL'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  upper_medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=100
  weight_prior=0
  plot_frequency=1

if calc_id == 'PLOT_ERF_1850_2005_2015':
  prior_id='UL'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  upper_medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=100
  weight_prior=0
  plot_frequency=1
  period_1850_2005_2015=True

if calc_id == 'PLOT_ERF_1750_2018':
  prior_id='UL'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  upper_medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=100
  weight_prior=0
  plot_frequency=1
  period_1750_2018=True

if calc_id == 'PLOT_ERF_AR5':
  # Use ='AR5_extended' aerosol forcing
  hist_forcing_version='AR5_extended'
  prior_id='UL'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  upper_medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=100
  weight_prior=0
  plot_frequency=1

if calc_id == 'PLOT_ERF_AR5_1750_2018':
  # Use ='AR5_extended' aerosol forcing
  hist_forcing_version='AR5_extended'
  prior_id='UL'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  upper_medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=100
  weight_prior=0
  plot_frequency=1
  period_1750_2018=True

if calc_id == 'PLOT_ERF_BELLCON':
  # Use Bellouin constrained aerosol forcing
  hist_forcing_version='Bellouin_2020_constrained'
  prior_id='UL'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  upper_medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=100
  weight_prior=0
  plot_frequency=1

if prior_id == 'US1':
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  pattern_effect_dependence2=0
  small_sample=1
  tiny_sample=0
  ntop=10000
  weight_prior=0
  plot_frequency=ntop/10

if prior_id == 'ULI':
  restrict_li_range=1
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  small_sample=1
  tiny_sample=0
  ntop=50000
  weight_prior=0
  plot_frequency=100

if calc_id == 'ULI_MEDIUM_SAMPLE':
  prior_id = 'ULI'
  restrict_li_range=1
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  small_sample=0
  medium_sample=1
  tiny_sample=0
  ntop=20000
  weight_prior=0
  plot_frequency=100

if calc_id == 'ULI_MEDIUM_SAMPLE_TECS':
  calc_tecs=True
  prior_id = 'ULI'
  restrict_li_range=1
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  small_sample=0
  medium_sample=1
  tiny_sample=0
  ntop=20000
  weight_prior=0
  plot_frequency=100

if calc_id == 'ULI_MEDIUM_SAMPLE_FAT':
  fat_tails=True
  prior_id = 'ULI'
  restrict_li_range=1
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  small_sample=0
  medium_sample=1
  tiny_sample=0
  ntop=20000
  weight_prior=0
  plot_frequency=100

if calc_id == 'US2_FULL':
  prior_id='US2'
  # alternative calculation without restricting li range and using
  # rejection sampling to support correlations.
  # same as US2 but restrict_li_range=False and weight_prior=0
  # same as US2_UNRESTRICTED_LI but weight_prior=0

  restrict_li_range=False
  weight_prior=0

  input_prior='ULI'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  tiny_sample=0
  ntop=1 # attempted
  ntop=50000 # attempted
  plot_frequency=2
  small_sample=1
  medium_sample=0
  upper_medium_sample=0

if calc_id == 'US2_FULL_PRIOR':
  prior_id='US2'
  # alternative calculation without restricting li range and using
  # rejection sampling to support correlations.
  restrict_li_range=False
  input_prior='ULI'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  tiny_sample=0
  ntop=1 # attempted
  ntop=50000 # attempted
  weight_prior=0
  plot_frequency=10
  small_sample=1
  medium_sample=0
  upper_medium_sample=0
  no_hist=1
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1

if calc_id =='US2':
  prior_id='US2'
  restrict_li_range=1
  input_prior='ULI'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  medium_sample=0
  small_sample=1
  tiny_sample=0
  ntop=50000
  weight_prior=1
  plot_frequency=100

if calc_id =='US2_MEDIUM_SAMPLE':
  prior_id='US2'
  restrict_li_range=1
  input_prior='ULI'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=20000
  weight_prior=1
  plot_frequency=100

if calc_id =='US2_MEDIUM_SAMPLE_TECS':
  calc_tecs=True
  prior_id='US2'
  restrict_li_range=1
  input_prior='ULI'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=20000
  weight_prior=1
  plot_frequency=100

if calc_id =='US2_MEDIUM_SAMPLE_FAT':
  fat_tails=True
  prior_id='US2'
  restrict_li_range=1
  input_prior='ULI'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=20000
  weight_prior=1
  plot_frequency=100

if calc_id =='US2_PROCBU':
  no_hist=1
  no_paleo_cold=1
  no_paleo_hot=1
  no_hist=1
  no_process_bu=0
  no_process_ec=1
  plot_process_likelihood=1
  prior_id='US2'
  restrict_li_range=1
  input_prior='ULI'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=50000
  weight_prior=1
  plot_frequency=100

if calc_id =='US2_UNRESTRICTED_LI':
  restrict_li_range=False
  weight_prior=1

  prior_id='US2'
  input_prior='ULI'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  medium_sample=1
  small_sample=0
  tiny_sample=0
  ntop=50000
  plot_frequency=ntop/100
  plot_frequency=100

if calc_id =='US2_PRIOR':
  prior_id='US2'
  restrict_li_range=0
  input_prior='ULI'
  plot_ecs_pdf=1
  plot_priors=0
  quick_and_dirty=0
  pattern_effect_dependence=0
  medium_sample=0
  small_sample=1
  tiny_sample=0
  ntop=50000
  weight_prior=1
  plot_frequency=ntop/100
  no_hist=1
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0

if calc_id == 'PRIORS':
  plot_ecs_pdf=0
  plot_priors=1
  quick_and_dirty=0
  pattern_effect_dependence=0
  plot_li_likelihoods=1
  medium_sample=0
  prior_sample=1
  small_sample=0
  tiny_sample=0
  ntop=1
  restrict_li_range=0
  weight_prior=True
  plot_frequency=ntop/10

  #restricted sampling needs small ntop e.g. 10

if calc_id == 'CLASSIC_HIST':
  classic_hist_pdf=True
  no_pattern_effect=True
  # Prior is not actually used for this calculation but it's
  # easier to calculate than to modify the code to not do so
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  ntop=20
  upper_medium_sample=1
  plot_priors=0
  plot_ecs_pdf=1

if calc_id == 'CLASSIC_HIST_ALT1':
  # Table 4.1 2nd row
  classic_hist_pdf=True
  no_pattern_effect=True
  # Prior is not actually used for this calculation but it's
  # easier to calculate than to modify the code to not do so
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  ntop=20
  upper_medium_sample=1
  plot_priors=0
  plot_ecs_pdf=1
  alternative_hist_dt=True
  alternative_DT_o=0.96  # update 19/3/20
  alternative_DT_sd=0.14/1.64

if calc_id == 'CLASSIC_HIST_ALT2':
  # Table 4.1 3rd row
  classic_hist_pdf=True
  no_pattern_effect=True
  # Prior is not actually used for this calculation but it's
  # easier to calculate than to modify the code to not do so
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  ntop=20
  upper_medium_sample=1
  plot_priors=0
  plot_ecs_pdf=1
  alternative_hist_dt=True
  alternative_DT_o=1.02 # Updated 19/3/20
  alternative_DT_sd=0.12/1.64
  alternative_pi_period=True

if calc_id == 'CLASSIC_HIST_ALT3':
  # Use ='AR5_extended' aerosol forcing
  hist_forcing_version='AR5_extended'
  classic_hist_pdf=True
  no_pattern_effect=True
  # Prior is not actually used for this calculation but it's
  # easier to calculate than to modify the code to not do so
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  ntop=20
  upper_medium_sample=1
  plot_priors=0
  plot_ecs_pdf=1

if calc_id == 'CLASSIC_HIST_ALT4':
  # Use Bellouin constrained aerosol forcing
  hist_forcing_version='Bellouin_2020_constrained'
  classic_hist_pdf=True
  no_pattern_effect=True
  # Prior is not actually used for this calculation but it's
  # easier to calculate than to modify the code to not do so
  prior_id='US1'
  no_hist=0
  no_paleo_cold=1
  no_paleo_hot=1
  no_process_bu=1
  no_process_ec=1
  plot_process_likelihood=0
  ntop=20
  upper_medium_sample=1
  plot_priors=0
  plot_ecs_pdf=1

def normalise_likelihood(bin_centres,likelihood):
  tmp=np.copy(likelihood)
  #tmp[bin_centres > 10]=0
  #tmp[bin_centres < -10]=0
  normalised_likelihood=likelihood/tmp.max()
  return(normalised_likelihood)

def restrict_lambda_i_samples(s,l1,l2,l3,l4,l5,l6,f,prior_weights,nsd):

  (l1_mu,l1_sd,l2_mu,l2_sd,l3_mu,l3_sd,l4_mu,l4_sd,l5_mu,l5_sd,l6_mu,l6_sd)=process_lambda_i_parameters()

  keep=\
  (l1 > l1_mu-nsd*l1_sd) & (l1 < l1_mu+nsd*l1_sd) & \
  (l2 > l2_mu-nsd*l2_sd) & (l2 < l2_mu+nsd*l2_sd) & \
  (l3 > l3_mu-nsd*l3_sd) & (l3 < l3_mu+nsd*l3_sd) & \
  (l4 > l4_mu-nsd*l4_sd) & (l4 < l4_mu+nsd*l4_sd) & \
  (l5 > l5_mu-nsd*l5_sd) & (l5 < l5_mu+nsd*l5_sd) & \
  (l6 > l6_mu-nsd*l6_sd) & (l6 < l6_mu+nsd*l6_sd)

  s=s[keep]
  l1=l1[keep]
  l2=l2[keep]
  l3=l3[keep]
  l4=l4[keep]
  l5=l5[keep]
  l6=l6[keep]
  f=f[keep]
  prior_weights=prior_weights[keep]

  return(s,l1,l2,l3,l4,l5,l6,prior_weights,f)

def calc_uniform_s_prior_from_li(n_samples,normal01_2co2,weights,restrict_li_range,input_prior='ULI',weight_prior=False,full_importance_density_scaling=np.array([0])):

  # Start with uniform priors U(-10,10) on L1 - L5
  # Need to sample more points as many get thrown away below

  fin=f2co2_process(normal01_2co2)

  label='US2'

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  if full_importance_density_scaling.shape[0] != 1:

    importance_density_scaling=full_importance_density_scaling

  else:

    # Calculate importance_sensity_scaling
    # importance_sensity_scaling is an array containing a scaling factor for each bin
    # multiplying the ULI prior on S by this should give a uniform prior on S

    # sample from ULI prior

    boost_sample_factor=10

    for iter in range(0,1000):

      min=-10
      max=10
      l1=np.random.uniform(min,max,n_samples*boost_sample_factor)
      l2=np.random.uniform(min,max,n_samples*boost_sample_factor)
      l3=np.random.uniform(min,max,n_samples*boost_sample_factor)
      l4=np.random.uniform(min,max,n_samples*boost_sample_factor)
      l5=np.random.uniform(min,max,n_samples*boost_sample_factor)
      l6=np.random.uniform(min,max,n_samples*boost_sample_factor)

      normal02_2co2=np.random.normal(0,1,n_samples*boost_sample_factor)
      f=f2co2_process(normal02_2co2)

      l=l1+l2+l3+l4+l5+l6
      s=-1*f/l

      # make prior density for proposal distribution

      this_prior, bin_boundaries = np.histogram(s, bins=bin_boundaries)

      if iter == 0:
        prior=np.copy(this_prior)
      else:
        prior=prior+this_prior

      print ("creating importance density scaling iter=",iter)

      shape = prior.shape

    # calculate importance density scaling
    # the importance density scaling is a set of weights, one for each bin
    # in the histogram / pdf
    # scaling the prior by these weights should result in a uniform prior on S

    #      print ('importance_density_scaling.shape:')
    #      print (importance_density_scaling.shape)
    #      print ('bin_boundaries.shape:')
    #      print (bin_boundaries.shape)

    importance_density_scaling=np.copy(prior)
    importance_density_scaling[bin_centres < 0.0] = 1.0
    importance_density_scaling[bin_centres >= 20.0] = 1.0
    # limit scaling for 0 or very small cases
    # OK because resulting distribution will look non-uniform if there's a problem
    importance_density_scaling[importance_density_scaling == 0.0] = 1.0
    importance_density_scaling=1.0/importance_density_scaling
    importance_density_scaling[bin_centres < 0.0] = 0.0
    importance_density_scaling[bin_centres >= 20.0] = 0.0

    # normalise scaling to have max 1 to support rejection sampling
    # print ('np.max(importance_density_scaling)=',iter, np.max(importance_density_scaling))

    importance_density_scaling=importance_density_scaling/np.max(importance_density_scaling)


  # now sample from the ULI prior, throwing away any samples with s outside of [-20,20]
  # to improve sample efficiency.  This won't affect the prior or posterior because any
  # samples outside the range [0,20] would have their weights scaled by zero anyway.

  # if weight_prior=0, use rejection sampling to subsample to give a uniform prior on S
  # if weight_prior=1, calculate a weight for each sample by looking up the scaling for
  # its S value

  for iter in range(0,1000000):

    #f=np.copy(fin)
    normal02_2co2=np.random.normal(0,1,n_samples)
    f=f2co2_process(normal02_2co2)

    if input_prior == 'ULI':

      if restrict_li_range:

        # if we are restricting the li range we sample uniformly over a +/- 6 sd range
        # of the process likelihoods
        # This increases sample efficiency by ignoring values which will be weighted to ~zero by process BU likelihood
        # and hence don't contribute to the posterior

        nsd=6
        (l1_mu,l1_sd,l2_mu,l2_sd,l3_mu,l3_sd,l4_mu,l4_sd,l5_mu,l5_sd,l6_mu,l6_sd)=process_lambda_i_parameters()

        l1=np.random.uniform(l1_mu-nsd*l1_sd,l1_mu+nsd*l1_sd,n_samples)
        l2=np.random.uniform(l2_mu-nsd*l2_sd,l2_mu+nsd*l2_sd,n_samples)
        l3=np.random.uniform(l3_mu-nsd*l3_sd,l3_mu+nsd*l3_sd,n_samples)
        l4=np.random.uniform(l4_mu-nsd*l4_sd,l4_mu+nsd*l4_sd,n_samples)
        l5=np.random.uniform(l5_mu-nsd*l5_sd,l5_mu+nsd*l5_sd,n_samples)
        l6=np.random.uniform(l6_mu-nsd*l6_sd,l6_mu+nsd*l6_sd,n_samples)

      else:

        # otherwise sample from the ULI prior

        min=-10
        max=10
        l1=np.random.uniform(min,max,n_samples)
        l2=np.random.uniform(min,max,n_samples)
        l3=np.random.uniform(min,max,n_samples)
        l4=np.random.uniform(min,max,n_samples)
        l5=np.random.uniform(min,max,n_samples)
        l6=np.random.uniform(min,max,n_samples)

      l=l1+l2+l3+l4+l5+l6
      s=-1*f/l

    # throw away values outside of the uniform S range [-20,20]
    # This makes the sampling more efficient but doesn't affect the posterior because these have a prior probability of zero

    discard_outside_s_plus_minus_20=True
    if discard_outside_s_plus_minus_20:

      smax=20
      smin=-20

      l1=l1[s >= smin]
      l2=l2[s >= smin]
      l3=l3[s >= smin]
      l4=l4[s >= smin]
      l5=l5[s >= smin]
      l6=l6[s >= smin]

      f=f[s >= smin]
      s=s[s >= smin]

      l1=l1[s < smax]
      l2=l2[s < smax]
      l3=l3[s < smax]
      l4=l4[s < smax]
      l5=l5[s < smax]
      l6=l6[s < smax]

      f=f[s < smax]
      s=s[s < smax]

    n_sampled=s.shape[0]

    #for each sampled s, calculate the index sbins of the bin it falls into

    sbins = ( s - bin_boundaries[0] ) / bin_width
    sbins = sbins.astype(int)

    #for each sampled s, calculate the scaling using the value for the relevant bin
    scalings=importance_density_scaling[sbins]

    if weight_prior == 0:

      #print ('before rejection sampling code:')
      #print ('iter=',iter)
      #print ('s.shape=',s.shape)

      #sample a uniform value u [0,1]

      u=np.random.uniform(0,1,n_sampled)

      # if u < importance_density_scaling(ibin) then accept the sample else reject it

      #      print ('s.shape=',s.shape)
      #      print ('scalings.shape=',scalings.shape)
      #      print ('u.shape=',u.shape)
      #      print ('sbins.shape=',sbins.shape)
      #
      #      print ('scalings=',scalings)
      #      print ('u=',u)
      #      print ('sbins=',sbins)

      # print (s.shape,u.shape,scalings.shape)
      s=s[u <= scalings]
      l1=l1[u <= scalings]
      l2=l2[u <= scalings]
      l3=l3[u <= scalings]
      l4=l4[u <= scalings]
      l5=l5[u <= scalings]
      l6=l6[u <= scalings]
      f=f[u <= scalings]

      # now we have a sample that's uniform in S

      # optionally throw out samples which will contribute little
      # weight to posterior because of 6D process likelihood.
      # it turns out this is very inefficient for ULI/US2 so will
      # have to do something else in this case.

    #print ('before call to restrict_lambda_i_samples')
    #print ('iter=',iter)
    #print ('s.shape=',s.shape)

    if (iter == 0):
      sall=np.copy(s)
      l1all=np.copy(l1)
      l2all=np.copy(l2)
      l3all=np.copy(l3)
      l4all=np.copy(l4)
      l5all=np.copy(l5)
      l6all=np.copy(l6)
      fall=np.copy(f)
    else:
      sall=np.append(sall,s)
      l1all=np.append(l1all,l1)
      l2all=np.append(l2all,l2)
      l3all=np.append(l3all,l3)
      l4all=np.append(l4all,l4)
      l5all=np.append(l5all,l5)
      l6all=np.append(l6all,l6)
      fall=np.append(fall,f)

    (datetime.datetime.now())
    #if iter % 10 == 1:
    if 1:
      print ('iter=',iter,'s.shape[0]',s.shape[0],'sall.shape[1]',sall.shape[0],'n_samples=',n_samples,'progress =',sall.shape[0]*100.0/n_samples,'%')
      print (datetime.datetime.now())
    if sall.shape[0] >= n_samples:
      print ('iter=',iter,'sall.shape (1)',sall.shape[0],'n_samples=',n_samples,'progress =',sall.shape[0]*100.0/n_samples,'%')
      print ('calc_uniform_s_prior_from_li: break!')
      break

  print ('loop finished: sall.shape[0]=',sall.shape[0],l1all.shape[0],
        l1all.shape[0],l2all.shape[0],l3all.shape[0],
        l4all.shape[0], l5all.shape[0], l6all.shape[0],fall.shape[0])

  # reduce to required number of samples
  f=fall[0:n_samples]
  s=sall[0:n_samples]
  l1=l1all[0:n_samples]
  l2=l2all[0:n_samples]
  l3=l3all[0:n_samples]
  l4=l4all[0:n_samples]
  l5=l5all[0:n_samples]
  l6=l6all[0:n_samples]

  prior, bin_boundaries = np.histogram(s, bins=bin_boundaries)
  shape = prior.shape

  if weight_prior:

    n_sampled=s.shape[0]

    #for each sampled s, calculate the index sbins of the bin it falls into

    sbins = ( s - bin_boundaries[0] ) / bin_width
    sbins = sbins.astype(int)

    #for each sampled s, calculate the scaling using the value for the relevant bin

    scalings=importance_density_scaling[sbins]

    # use scalings as importance weights for prior
    prior_weights=np.copy(scalings)

  else:

    # set equal prior weights if we used rejection sampling

    prior_weights=np.full(n_samples,1.0)

  #prior=normalise_pdf(prior,bin_width)

  return(prior,f,l1,l2,l3,l4,l5,l6,s,prior_weights,importance_density_scaling)

def ecs_cdf(count_in,bin_boundaries,bin_width,label):

  count=normalise_pdf(count_in,bin_width)

  if label == 'Loopsum Posterior':
    kernel_smoothing=False

    if kernel_smoothing:

      # use kernel smoothing for S

      pdf=normalise_pdf(count_in,bin_width)

      # apply Gausian Kernel smoothing to pdf
      # JDA uses sd of 0.1
      kernel_sd = 0.1

      smoothed_pdf=np.copy(pdf)

      # apply kernel filter over range 0-20
      istart=int(n_bins/2)
      iend=int(n_bins/2+2000)

      for i in range(istart,iend):
        x=bin_centres[i]
        k = np.exp(-1*( x - bin_centres ) ** 2 / (2 * kernel_sd ** 2))
        smoothed_pdf[i] = np.sum(pdf * k)
        #print ('i=',i,'bin_center=',bin_centres[i],x,posterior[i],smoothed_posterior[i])

      count=normalise_pdf(smoothed_pdf,bin_width)

    else:
      count=normalise_pdf(count_in,bin_width)

  slabel=label
  cdf=np.float64(0.0)
  ecs_mean=np.float64(0.0)
  ecs_mode=np.float64(0.0)
  ecs_weight=np.float64(0.0)
  # ubin_boundaries contains upper limits of bin_boundaries
  ubin_boundaries=bin_boundaries+bin_width
  shape = count.shape
  ecs5=0.0
  ecs10=0.0
  ecs17=0.0
  ecs20=0.0
  ecs25=0.0
  ecs50=0.0
  ecs75=0.0
  ecs80=0.0
  ecs83=0.0
  ecs90=0.0
  ecs95=0.0
  p_ecs_gt_3=0.0

  for x in range(0, shape[0]):
           bin_centre=ubin_boundaries[x]+bin_width/2
           cdf=cdf+count[x]*bin_width*0.5
           ecs_mean=ecs_mean+bin_centre*count[x]
           ecs_weight=ecs_weight+count[x]
           if bin_centre < 2:
                 p_ecs_lt_2=cdf*100
           if bin_centre < 1:
                 p_ecs_lt_1=cdf*100
           if bin_centre < 1.5:
                 p_ecs_lt_1p5=cdf*100
           if bin_centre < 0:
                 p_ecs_lt_0=cdf*100
           if bin_centre > 3 and p_ecs_gt_3 == 0:
             p_ecs_gt_3=cdf*100
           if bin_centre < 3:
             p_ecs_lt_3=cdf*100
           if bin_centre < 4:
             p_ecs_lt_4=cdf*100
           if bin_centre < 4.5:
             p_ecs_lt_4p5=cdf*100
           if bin_centre < 6.0:
             p_ecs_lt_6=cdf*100
           if cdf > 0.05 and ecs5 == 0:
             ecs5=bin_boundaries[x]+bin_width/2
           if cdf > 0.17 and ecs17 == 0:
             ecs17=bin_boundaries[x]+bin_width/2
           if cdf > 0.25 and ecs25 == 0:
             ecs25=bin_boundaries[x]+bin_width/2
           if cdf > 0.1 and ecs10 == 0:
             ecs10=bin_boundaries[x]+bin_width/2
           if cdf > 0.2 and ecs20 == 0:
             ecs20=bin_boundaries[x]+bin_width/2
           if cdf > 0.5 and ecs50 == 0:
             ecs50=bin_boundaries[x]+bin_width/2
           if cdf > 0.75 and ecs75 == 0:
             ecs75=bin_boundaries[x]+bin_width/2
           if cdf > 0.83 and ecs83 == 0:
             ecs83=bin_boundaries[x]+bin_width/2
           if cdf > 0.8 and ecs80 == 0:
             ecs80=bin_boundaries[x]+bin_width/2
           if cdf > 0.90 and ecs90 == 0:
             ecs90=bin_boundaries[x]+bin_width/2
           if cdf > 0.95 and ecs95 == 0:
             ecs95=bin_boundaries[x]+bin_width/2
           cdf=cdf+count[x]*bin_width*0.5
  p_ecs_gt_4p5=100-p_ecs_lt_4p5

  n=np.shape(bin_boundaries)
  n=n[0]-1
  ecs_mean=ecs_mean/ecs_weight

  ecs_mode_index=np.argmax(count)
  ecs_mode=bin_boundaries[ecs_mode_index]+bin_width/2

  print ('bin_width = ',bin_width)

  print ( 'Prior;P(ECS <1);P(ECS <1.5);P(ECS <2);P(ECS >4);P(ECS >4.5);P(ECS >6);P(1.5 -4.5);5th pile;10th pile;17th pile;20th pile;25th pile;50th pile;75th pile;80th pile;83rd pile;90th pile;95th pile;Mode;Mean' )
  qlabel=label
  f='10.2f'
  f=''
  print ("qlabel=",qlabel, p_ecs_lt_1, p_ecs_lt_1p5, p_ecs_lt_2,100-p_ecs_lt_4,100-p_ecs_lt_4p5, 100-p_ecs_lt_6,p_ecs_lt_4p5-p_ecs_lt_1p5,ecs5,ecs10, ecs17, ecs20,ecs25,ecs50, ecs75,ecs80,ecs83, ecs90, ecs95,ecs_mode,ecs_mean)

  print

  #  print ('P(ECS < 0) = ' , p_ecs_lt_0
  #  print ('P(ECS < 1.5) = ' , p_ecs_lt_1p5
  #  print ('P(ECS < 2.0) = ' , p_ecs_lt_2
  #  print ('P(ECS < 3.0) = ' , p_ecs_lt_3
  #  print ('P(ECS > 3.0) = ' , p_ecs_gt_3
  #  print ('P(ECS > 4.0) = ' , 1-p_ecs_lt_4
  #  print ('P(ECS > 4.5) = ' , 1-p_ecs_lt_4p5
  #  print ('P(ECS > 6.0) = ' , 1-p_ecs_lt_6
  #  print ('P(1.5 < ECS < 4.5) = ' , p_ecs_lt_4p5 - p_ecs_lt_1p5
  #  print ('5th percentile = ' , ecs5
  #  print ('50th percentile = ' , ecs50
  #  print ('95th percentile = ' , ecs95
  #  print ('5th minus 50th = ' , ecs5-ecs50
  #  print ('95th minus 50th = ' , ecs95-ecs50
  #print ('ptest' , ptest
  print ('histogram total' , cdf)
  print ('cdf=', cdf)
  assert (cdf > 0.99)
  assert (cdf <= 1.01)

  print (format(qlabel,'50s'),format('5-95% range:','10s'),format(ecs5 ,f),'-',format(ecs95 ,f))

  slabel2=slabel+'.txt'
  slabel2=slabel2.replace(' ','_')
  slabel2=slabel2.replace('-','_')
  slabel2=slabel2.replace('/','_')

  #file = open(slabel2,"w")
  #file.write(bin_boundaries)
  #file.write(count)
  #file.close()

  return (p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode)

def read_hist_forc():

  #start by making the forcing according to Piers' code and commentary

  #note the only edit I made to his files is to make the col names automatically readable for the 2nd one (ie removing spaces)

  # Piers Forster's code retained in comments marked PF
  # James Annan's code retained in comments marked JDA

  import numpy as np
  import pandas as pd

  # Read in data

  from io import StringIO   # StringIO behaves like a file object
  #from google.colab import drive
  #drive.mount('/content/gdrive')

  # **** Read in historical forcing upper, mid and lower ranges ****

  if hist_forcing_version == 'Bellouin_2020_eqn_8':
    hist_forcing_ranges='WCRP_forcing_ranges_200311.csv'
    hist_monthly_forcings='Forster_forcings_2018_WCRP_Bellouinaerosol2.csv'

  if hist_forcing_version == 'Bellouin_2020_constrained':
    hist_forcing_ranges='WCRP_forcing_ranges_Bellouinaerosol_constrained_200321.csv'
    hist_monthly_forcings='Forster_forcings_2018_WCRP_Bellouinaerosol2.csv'

  if hist_forcing_version == 'AR5_extended':
    hist_forcing_ranges='IPCCAR5_forcing_ranges_200320.csv'
    hist_monthly_forcings='Forster_forcings_2018_WCRP_AR5aerosol_200320.csv'

  best_upper_lower = pd.read_csv(inpath+hist_forcing_ranges,nrows=3)

  best=best_upper_lower.iloc[0]
  lower=best_upper_lower.iloc[1]
  upper=best_upper_lower.iloc[2]

  # PF ; assumes 90% are 1.65 sigma
  # PF sdlower=(best-lower)/1.65
  # PF sdupper=(upper-best)/1.65

  # MJW calulate standard deviations which would give lower and upper values as 5-95% ranges
  # assumes 90% are 1.64 sigma
  sdlower=(best-lower)/1.64
  sdupper=(upper-best)/1.64

  print ("best:")
  print (best)
  print ("lower:")
  print (lower)
  print ("upper:")
  print (upper)

  print ("sdlower:")
  print (sdlower)
  print ("sdupper:")
  print (sdupper)

  # **** Read in AR5 annual mean forcing data ****

  df=pd.read_csv(inpath+hist_monthly_forcings,delimiter=',',skiprows=19)

  #print (forc_table["Year/month"]
  df = df.set_index(['Year/month'])

  df_baseline=df.loc[1861:1880]
  #df_recent=df.loc[2002:2017]
  df_recent=df.loc[2006:2018] # updated 16/3/20

  if alternative_pi_period:
    df_baseline=df.loc[1850:1900]

  if period_1850_2005_2015:
    df_baseline=df.loc[[1850,1850]]
    df_recent=df.loc[2005:2015]

  if period_1750_2018:
    df_baseline=df.loc[[1750,1750]]
    df_recent=df.loc[[2018,2018]]

  print (df_baseline)
  print (df_recent)

  # JDA: CO2 forcing is based on F2x so I will back-calculate the doubling fraction from this
  #
  # JDA: now generate a sample of all the other forcings

  uncs_neg = ( best - lower ) / (1.64 * best )
  uncs_pos = ( upper - best ) / (1.64 * best )

  return (df_recent,df_baseline,uncs_neg,uncs_pos)

def hist_forc(nsamp,df_recent,df_baseline,uncs_neg,uncs_pos):

  f_other =np.zeros(nsamp)

  # JDA: incs = 2:12
  # JDA: incs = incs[-4] #get rid of RFari
  # JDA: for (i in incs)
  # JDA: {
  # JDA: print(i)

  # PF for imc=0,999 do begin
  # PF ; choose 12 random fractions of each forcing term
  # MJW This comment would seem to be wrong - Piers's code uses 11 forcing
  # components because it (correctly) ignores RFari
  # PF p=randomn(seed,12)
  # PF frac=p
  # PF ilow=where (p le 0.0)
  # PF ihigh=where (p ge 0.0)
  # PF frac[ilow]=(best[ilow]+p[ilow]*sdlower[ilow])/best[ilow]
  # PF frac[ihigh]=(best[ihigh]+p[ihigh]*sdupper[ihigh])/best[ihigh]

  # PF ar5total[imc,*]=ar5_85[0,*]*frac[0]+ar5_85[1,*]*frac[1]+ar5_85[2,*]*frac[2]
  #                   +ar5_85[3,*]*frac[3]+ar5_85[5,*]*frac[5]+ar5_85[6,*]*frac[6]
  #                   +ar5_85[7,*]*frac[7]+ar5_85[8,*]*frac[8]+ar5_85[9,*]*frac[9]
  #                   +ar5_85[10,*]*frac[10]+ar5_85[11,*]*frac[11]

  labels = ['Other WMGHG','O3 (T)','O3(S)','total aerosol ERF','ERF LUC',
            'Vapour', 'BC snow','contrails','Solar','Volcanic']

  # loop over non-co2 (other) forcings
  # NOTE we don't include RFari column in data files
  # because this is included in the total aerosol ERF

  for label in labels:

    # JDA: for each forcing, generate a vector of deviates and scale according to uncertainties

    # JDA devs = rnorm(nsamp)

    baseline=df_baseline[label].to_numpy()
    baseline=np.mean(baseline)

    recent=df_recent[label].to_numpy()
    recent=np.mean(recent)

    print (label,"baseline:",baseline)
    print (label,"recent:",recent)
    print (label,"recent-baseline:",recent-baseline)

    f_bin_boundaries=np.linspace(-10,10,n_bins+1)
    f_bin_width=f_bin_boundaries[1]-f_bin_boundaries[0]
    f_bin_centres=f_bin_boundaries[0:n_bins]+f_bin_width/2

    if hist_forcing_version == 'Bellouin_2020_eqn_8' and label == 'total aerosol ERF':
      # Calculate Bellouin et al 2020 eqn 8 aerosol forcing PDF instead of ranges file
      ERFaer = bellouin_ringberg_stats_fix_200311_ERFaer(nsamp)
      #devscaled=ERFaer/np.median(ERFaer)-1.0
      # rescale by the factor that reproduces Bellouin et al 2020 PDF between 1850 and 2005-2015
      peak_ERFaer = -0.87
      devscaled = ERFaer / peak_ERFaer - 1.0

      print ("np.median(ERFaer)=",np.median(ERFaer))

      flabel="ERFaer"
      f=np.copy(ERFaer)
      f_pdf, f_bin_boundaries = np.histogram(f, bins=f_bin_boundaries)
      ([p_f_lt_1p5,p_f_gt_4p5,f5,f95,f_mean,f_mode]) = ecs_cdf(f_pdf,f_bin_boundaries,f_bin_width,flabel)
      print(flabel,"=",f_mean,f5,f95)

      flabel="devscaled"
      f=np.copy(devscaled)
      f_pdf, f_bin_boundaries = np.histogram(f, bins=f_bin_boundaries)
      ([p_f_lt_1p5,p_f_gt_4p5,f5,f95,f_mean,f_mode]) = ecs_cdf(f_pdf,f_bin_boundaries,f_bin_width,flabel)
      print(flabel,"=",f_mean,f5,f95)

    else:

      devs = np.random.normal(loc=0.0, scale=1.0, size=nsamp)

      print (label,"uncs_pos:",uncs_pos[label])
      print (label,"uncs_neg:",uncs_neg[label])

      # JDA: devscaled = ((devs > 0)*devs*as.matrix(uncs_pos)[i]) + ((devs < 0)*devs*as.matrix(uncs_neg)[i])

      devscaled = (devs >= 0) * devs * uncs_pos[label]
      devscaled = devscaled + (devs < 0) * devs * uncs_neg[label]

    ## JDA: then use these deviates as a _scaling_ on the actual forcing change

    # JDA: f_other = f_other + (mean(forc_table[recent,i+1]) - mean(forc_table[baseline,i+1])) * (1+devscaled)

    print ( "devscaled.shape",devscaled.shape)

    f_component = (recent - baseline) * (1+devscaled)
    f_other = f_other + f_component

    if label == 'total aerosol ERF':
      total_aerosol_erf=np.copy(f_component)

    # write out stats for forcing components

    f=np.copy(f_component)
    flabel=label
    f_pdf, f_bin_boundaries = np.histogram(f, bins=f_bin_boundaries)
    ([p_f_lt_1p5,p_f_gt_4p5,f5,f95,f_mean,f_mode])=ecs_cdf(f_pdf,f_bin_boundaries,f_bin_width,flabel)
    print(flabel,"=",f_mean,f5,f95)

  # JDA: CO2_doub = (mean(forc_table[recent,2]) - mean(forc_table[baseline,2])) / np.mean(F2x)

  # MJW: I do it differently from James.  I just return the mean hist CO2 forcing.
  # This code doesn't then need to know what the 2CO2 forcing is.

  label='CO2'
  baseline=df_baseline[label].to_numpy()
  baseline=np.mean(baseline)

  recent=df_recent[label].to_numpy()
  recent=np.mean(recent)

  CO2_hist = (recent - baseline)

  print ("CO2_hist:",CO2_hist)

  print ("f_other mean:",np.mean(f_other))

  print ("f_other:",f_other)

  return (f_other,CO2_hist,total_aerosol_erf)

def process_likelihood_weight_l_pattern_effect(l,lam_p_p):

  (lambda_mu,lambda_sd)=process_lambda_parameters()

  lam_p_p_mu=np.mean(lam_p_p)
  lam_p_p_sd=np.std(lam_p_p)

  mu=lambda_mu+lam_p_p_mu
  sd=np.sqrt(lambda_sd**2+lam_p_p_sd**2)

  weight=dnorm(mu,l,sd)

  return (weight)

def process_likelihood_weight_l(l):

  (lambda_mu,lambda_sd)=process_lambda_parameters()

  if calc_tecs:
    # use transfer function if calculating true equilibrium climate sensitivity (TECS)
    # in this case prior l is lambda_TECS
    # TECS = S * (1+transfer)
    # lambda_TECS = lambda_S / (1+transfer)
    # lambda_S = lambda_TECS * (1+transfer)
    # process evidence is a function of lambda_S
    # so we want to use lambda_S = l*(1+transfer)
    if fat_tails:
      weight=dfat(lambda_mu,l*(1+transfer),lambda_sd)
    else:
      # baseline TECS
      weight=dnorm(lambda_mu,l*(1+transfer),lambda_sd)
  else:
    if fat_tails:
      weight=dfat(lambda_mu,l,lambda_sd)
    else:
      # baseline S
      weight=dnorm(lambda_mu,l,lambda_sd)

  return (weight)

#def lambda_process():
#  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()
#
#  femu,fesigma = process_lambda_parameters()
#  fes=np.random.normal(femu,fesigma,n_samples)
#
#  return fes

def process_lambda_parameters():

  lambda_mu=-1.30 # Updated 6/3/20
  lambda_sd= 0.44  # Updated 6/3/20

  return (lambda_mu,lambda_sd)

def process_lambda_i_parameters():

  l1_mu = -3.2 # planck
  l1_sd = 0.1
  l2_mu = 1.15 # LR+WV
  l2_sd = 0.15
  l3_mu = 0.3 # alb
  l3_sd = 0.15
  l4_mu = l_cld_mu # 0.45  # cloud Updated 6/3/20
  l4_sd = l_cld_sd # 0.33
  l5_mu = 0.   # Atmos Composition
  l5_sd = 0.15
  l6_mu = 0.   # Stratospheric feedback added 20/2/19
  l6_sd = 0.1

  return (l1_mu,l1_sd,l2_mu,l2_sd,l3_mu,l3_sd,l4_mu,l4_sd,l5_mu,l5_sd,l6_mu,l6_sd)

def process_likelihood_weight_li(l1,l2,l3,l4,l5,l6,plot_lis):

  #function to calculate the process 6D likelihood

  (l1_mu,l1_sd,l2_mu,l2_sd,l3_mu,l3_sd,l4_mu,l4_sd,l5_mu,l5_sd,l6_mu,l6_sd) \
                                            = process_lambda_i_parameters()

  print ("XX5 l5.shape=",l5.shape)

  if calc_tecs:
    # use transfer function if calculating true equilibrium climate sensitivity (TECS)
    # in this case priors l1-l6 are priors on lambda components of lambda_TECS
    # TECS = S * (1+transfer)
    # lambda_TECS = lambda_S / (1+transfer)
    # lambda_S = lambda_TECS * (1+transfer)
    # process evidence is a function of the six components of lambda_S
    # We need to make an assumption for how to distribute the scaling among
    # components.  We simply scale all components by the same fraction.
    # so we to use lambda_S_i = lambda_TECS_i*(1+transfer)
    if fat_tails:
      # fat tails
      like_1 = dfat(l1_mu,l1*(1+transfer),l1_sd)
      like_2 = dfat(l2_mu,l2*(1+transfer),l2_sd)
      like_3 = dfat(l3_mu,l3*(1+transfer),l3_sd)
      like_4 = dfat(l4_mu,l4*(1+transfer),l4_sd)
      like_5 = dfat(l5_mu,l5*(1+transfer),l5_sd)
      like_6 = dfat(l6_mu,l6*(1+transfer),l6_sd)
    else:
      # baseline tecs
      like_1 = dnorm(l1_mu,l1*(1+transfer),l1_sd)
      like_2 = dnorm(l2_mu,l2*(1+transfer),l2_sd)
      like_3 = dnorm(l3_mu,l3*(1+transfer),l3_sd)
      like_4 = dnorm(l4_mu,l4*(1+transfer),l4_sd)
      like_5 = dnorm(l5_mu,l5*(1+transfer),l5_sd)
      like_6 = dnorm(l6_mu,l6*(1+transfer),l6_sd)
  else:
    if fat_tails:
      # fat tails
      like_1 = dfat(l1_mu,l1,l1_sd)
      like_2 = dfat(l2_mu,l2,l2_sd)
      like_3 = dfat(l3_mu,l3,l3_sd)
      like_4 = dfat(l4_mu,l4,l4_sd)
      like_5 = dfat(l5_mu,l5,l5_sd)
      like_6 = dfat(l6_mu,l6,l6_sd)
    else:
      # baseline
      like_1 = dnorm(l1_mu,l1,l1_sd)
      like_2 = dnorm(l2_mu,l2,l2_sd)
      like_3 = dnorm(l3_mu,l3,l3_sd)
      like_4 = dnorm(l4_mu,l4,l4_sd)
      like_5 = dnorm(l5_mu,l5,l5_sd)
      like_6 = dnorm(l6_mu,l6,l6_sd)

  print ("like_1",like_1.shape)
  print ("like_2",like_2.shape)
  print ("like_3",like_3.shape)
  print ("like_4",like_4.shape)
  print ("like_5",like_5.shape)

  likelihood = like_1 * like_2 * like_3 * like_4 * like_5 * like_6

  if plot_lis and toploop_index % plot_frequency ==0:
    (lambda_mu,lambda_sd) = process_lambda_parameters()
    #plt.close('all')
    pos=2
    l=l1+l2+l3+l4+l5+l6

    plot_li('Planck_Feedback',l1_mu,l1,l1_sd,like_1,'r',1)
    plot_li('Lapse_rate_+_Water_Vapour_Feedback',l2_mu,l2,l2_sd,like_2,'g',1)
    plot_li('Surface_Albedo_Feedback',l3_mu,l3,l3_sd,like_3,'b',1)
    plot_li('Cloud_Feedback',l4_mu,l4,l4_sd,like_4,'c',1)
    plot_li('Atmospheric_Composition',l5_mu,l5,l5_sd,like_5,'m',1)
    plot_li('Stratospheric_Feedback',l6_mu,l6,l6_sd,like_6,'y',1)
    plot_li('Total_Lambda',lambda_mu,l,lambda_sd,likelihood,'k',6)

    plot_li_scat('l1_vs_l2',l1,l2)
    plot_li_scat('l1_vs_l2+l3+l4+l5+l6',l1,l2+l3+l4+l5+l6)

    plt.savefig(out_dir+'/'+calc_id+'1f.png')
    plt.close('all')
    print ('xx li likelihoods',toploop_index)

  return(likelihood)

def plot_li(label,mu,prior,sd,weights,color,pos):
  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()
  equal_weights=np.full(n_samples,1.0)
  likelihood_sample=np.random.normal(mu,sd,n_samples)
  unweighted_li_prior=sample2histogram(prior,bin_boundaries,equal_weights)
  li_prior=sample2histogram(prior,bin_boundaries,prior_weights)
  li_likelihood=sample2histogram(likelihood_sample,bin_boundaries,equal_weights)
  li_posterior=sample2histogram(prior,bin_boundaries,weights)

  # scale posterior to normalise to max in the plotting range
  tmp=np.copy(li_posterior)
  #tmp[bin_centres > 10]=0
  #tmp[bin_centres < -10]=0
  li_posterior=li_posterior/tmp.max()

  unweighted_li_prior_plot=normalise_pdf(unweighted_li_prior,bin_width)
  li_prior_plot=normalise_pdf(li_prior,bin_width)
  li_posterior_plot=normalise_pdf(li_posterior,bin_width)
  li_posterior_plot=normalise_pdf(li_posterior,bin_width)
  li_likelihood_plot=normalise_pdf(li_likelihood,bin_width)

  print ('plot_li label=',label,mu,sd)

  # **** plot l prior predictive distributions *****

  xpos=2

  plt.close('all')

  plt.xlim([-11,11])
  plt.ylim(0,2)
  if not plot_posterior:
    plt.ylim(0,1)

  color='k'
  plt.text(xpos,yrange[1]*1-0.1*(pos+1),label,color=color)
  #plt.plot(bin_centres,li_prior,color)
  #plt.plot(bin_centres,li_likelihood,color,linestyle='dashed')
  #plt.plot(bin_centres,li_posterior,color,linestyle='dotted')
  plt.plot(bin_centres,li_prior_plot,'g')
  plt.plot(bin_centres,unweighted_li_prior_plot,'r')
  plt.plot(bin_centres,li_likelihood_plot,'b')
  plt.plot(bin_centres,li_posterior_plot,'k')

  plt.savefig(out_dir+'/'+calc_id+'1_'+label+'.png')
  plt.close('all')

def plot_li_scat(label,x,y):

  xpos=2
  pos=7

  plt.close('all')

  plt.xlim([-11,11])
  plt.ylim([-11,11])

  plt.scatter(x,y)

  plt.text(xpos,yrange[1]*1-0.1*(pos+1),label)

  plt.savefig(out_dir+'/'+calc_id+'.scat.'+label+'.png')
  plt.close('all')

def plot_sample(sample,label):
  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()
  equal_weights=np.full(n_samples,1.0)
  pdf=sample2histogram(sample,bin_boundaries,equal_weights)

  plt.xlim([-5,5])
  plt.ylim(0,2)
  plt.title(label)
  plt.plot(bin_centres,pdf)
  print ('plot_sample',label)

  return (weight)

def process_likelihood_weight_emergent_constriants_s(l,f2x,lambda_ec=False):

  (mu_ec,sigma_ec) = emergent_constraints_lambda_parameters()

  weight=dnorm(mu_ec,l,sigma_ec)

  return (weight)

def paleo_cold_likelihood(l,F2x,transfer,F_prime,alpha):

  #function to calculate the paleo likelihood based on cold periods

  #change 30 Jan 2020 to account for 284ppm and Etminan = 0.57 of a doubling
  #unchanged in 16 March as F' comes from the synthesis code!

  #model for temperature change

  #DT = sum (0.57 F2x + F' ) / (l + alp DT)

  #note additive term not multiplicative for the alpha

  #where F' is the other (non-CO2) forcings

  #needs to be reorganised into a quadratic
  #note that I am integrating this equation to get the temp
  #as applied forcing df = int_0^dt of (l + alp dt) = l.dt + 1/2 alp.dt^2

  # Note l and lam have different sign conventions in this code
  # here l is negative and lam is positive
  # In James Annan's original R code input variable lambda and lam were
  # both positive
  # Hence -1 in this equation which isn't in original R code

  if calc_tecs:
    # don't use transfer function if calculating true equilibrium climate sensitivity
    lam = -1*l
  else:
    # lam is associated with quasi equilibrium
    # s_eqm = s * (1+transfer) so lam = - l / (1+transfer)
    lam = -1*l/(1+transfer)

  #using etminan et al forcing number (ratio)

  DT = (-1*lam + np.sqrt(lam**2 + 2*alpha*(F_prime - 0.57 * F2x)))/(alpha)

  DT_o = -5.
  DT_sd = 1.

  likelihood = dnorm(DT_o,DT,DT_sd)

  #fix any NaNs which are due to no roots = unbounded temp change
  likelihood[np.isnan(likelihood)] = 0.

  return(likelihood)

#making a paleo likelihood

#warm periods calculation

def paleo_hot_likelihood(l,F2x,transfer,ppm,ch4_fac,scalefac):

  #fixing at 284 now 30 Jan 2020
  #NB the ch4 factor already includes the F2x scaling so no change is needed,
  #this already implicitly accounts for any scaling of F2x.

  #DT = (ln(ppm/284)/ln(2) * (1+ch4_fac) * F2x * (1+scalefac) ) / l
  #where (1+scalefac) does the ESS/S ratio of about 1.4ish
  #ppm is the co2 eq level which may be 405 but probably a bit lower
  #

  DT_o = 3.
  DT_sd = 1.

  # Note l and lam have different sign conventions in this code
  # here l is negative and lam is positive
  # In James Annan's original R code input variable lambda and lam were
  # both positive
  # Hence -1 in this equation which isn't in original R code

  if calc_tecs:
    # don't use transfer function if calculating true equilibrium climate sensitivity (TECS)
    lam = -1*l
  else:
    # lam is associated with quasi equilibrium
    # s_eqm = s * (1+transfer) so lam = - l / (1+transfer)
    lam = -1*l/(1+transfer)

  DT = (np.log(ppm/284)/np.log(2) * (1+ch4_fac) * F2x * (1+scalefac) ) / lam

  likelihood = dnorm(DT_o,DT,DT_sd)

  #fix any NaNs which are due to no roots = unbounded temp change
  likelihood[np.isnan(likelihood)] = 0.

  return(likelihood)

def calc_prior_uniform_s(n_samples,normal01_2co2,weights):

  f=f2co2_process(normal01_2co2)

  label='s_uniform'

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  s=np.random.uniform(0,20,n_samples)

  l=-1*f/s

  prior=sample2histogram(s, bin_boundaries,weights)

  return(prior,f,l,s)

def sample2histogram(sin,bin_boundaries,weights,inspect_stats=None):

  #  from statsmodels.nonparametric.kde import KDEUnivariate

  s=np.copy(sin)
  kde=0
  n_bins=bin_boundaries.shape
  n_bins=n_bins[0]-1
  bin_width=bin_boundaries[1]-bin_boundaries[0]
  bin_centres=bin_boundaries[0:n_bins]+bin_width/2
  bin_lower_boundaries=bin_boundaries[0:n_bins-1]
  bin_upper_boundaries=bin_lower_boundaries+bin_width
  nsamp=s.shape
  nsamp=nsamp[0]
  if kde == 1:
    sub=100
    s=sin[0:sub]
    kde1= KDEUnivariate(sin)
    kde1.fit(weights=weights,
      bw=0.2,
      fft=False)
    pdf=[kde1.evaluate(xi) for xi in bin_centres]
  if kde == 2:
    kernel_half_width=0.2
    pdf=np.zeros(n_bins)
    for i in range(nsamp):
      kernel_index=((bin_lower_boundaries - kernel_half_width < s[i]).all() and (s[i] < bin_upper_boundaries + kernel_half_width).all())
      pdf[kernel_index]=pdf[kernel_index]+weights[i]
  if kde == 0:
    pdf, bin_boundaries = np.histogram(s, weights=weights, bins=bin_boundaries)

  #print('sample2histogram:',pdf)
  #pdf=normalise_pdf(pdf,bin_width)
  #print('sample2histogram:',pdf)

  if inspect_stats != None:
    # sort weights and select 20 largest values
    isort=(-weights).argsort()[:20]
    s_sorted=s[isort]
    weights_sorted=weights[isort]
    print("s_sorted",s_sorted)
    print("weights_sorted",weights_sorted)
    isort=(-pdf).argsort()[:20]
    pdf_sorted=pdf[isort]
    print("pdf_sorted",pdf_sorted)
    bin_centres_sorted=bin_centres[isort]
    print("bin_centres_sorted",bin_centres_sorted)

  return(pdf)

def calc_prior_uniform_lambda(n_samples,normal01_2co2,weights,lower=-20):

  f=f2co2_process(normal01_2co2)

  label='l_uniform'

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  l=np.random.uniform(lower,0,n_samples)

  s=-1*f/l

  prior=sample2histogram(s, bin_boundaries,weights)

  return(prior,f,l,s)

def calc_process_bu_prior(n_samples,normal01_2co2,weights):

  f=f2co2_process(normal01_2co2)

  label='l_normal'

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  (loc,scale)=process_lambda_parameters()

  l = np.random.normal(loc,scale,n_samples)

  s=-1*f/l

  prior=sample2histogram(s, bin_boundaries,weights)

  return(prior,f,l,s)

def calc_process_prior(n_samples,normal01_2co2,weights):

  f=f2co2_process(normal01_2co2)

  label='process_prior'

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  (lambda_mu,lambda_sd)=process_lambda_parameters()

  l=np.random.normal(lambda_mu,lambda_sd,n_samples)

  s=-1*f/l

  prior=sample2histogram(s, bin_boundaries,weights)

  return(prior,f,l,s)

def calc_prior_uniform_lambda_is(n_samples,normal01_2co2,weights,restrict_li_range):

  f = f2co2_process(normal01_2co2)

  label='li_uniform'

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  # if we are restricting the li range we sample uniformly over a +/- 4sd range
  # of the process likelihoods

  if restrict_li_range:

    nsd=6
    (l1_mu,l1_sd,l2_mu,l2_sd,l3_mu,l3_sd,l4_mu,l4_sd,
     l5_mu,l5_sd,l6_mu,l6_sd) = process_lambda_i_parameters()

    l1=np.random.uniform(l1_mu-nsd*l1_sd,l1_mu+nsd*l1_sd,n_samples)
    l2=np.random.uniform(l2_mu-nsd*l2_sd,l2_mu+nsd*l2_sd,n_samples)
    l3=np.random.uniform(l3_mu-nsd*l3_sd,l3_mu+nsd*l3_sd,n_samples)
    l4=np.random.uniform(l4_mu-nsd*l4_sd,l4_mu+nsd*l4_sd,n_samples)
    l5=np.random.uniform(l5_mu-nsd*l5_sd,l5_mu+nsd*l5_sd,n_samples)
    l6=np.random.uniform(l6_mu-nsd*l6_sd,l6_mu+nsd*l6_sd,n_samples)
  else:
    limit=10
    l1=np.random.uniform(-1*limit,limit,n_samples)
    l2=np.random.uniform(-1*limit,limit,n_samples)
    l3=np.random.uniform(-1*limit,limit,n_samples)
    l4=np.random.uniform(-1*limit,limit,n_samples)
    l5=np.random.uniform(-1*limit,limit,n_samples)
    l6=np.random.uniform(-1*limit,limit,n_samples)

  l = l1 + l2 + l3 + l4 + l5 + l6
  s = -1 * f / l
  prior = sample2histogram(s, bin_boundaries, weights)

  return(prior,f,l,s,l1,l2,l3,l4,l5,l6)

def setup_bins():

  # reducing this number biases the calculation of the probabilities
  n_samples=10000000
  n_bins=100000
  bin_boundaries=np.linspace(-1000,1000,n_bins+1)

  # smoothing experiment
  n_samples=10000000
  n_bins=1000
  bin_boundaries=np.linspace(-100,100,n_bins+1)

  # smoothing experiment 2
  #n_samples=500000
  #n_bins=1000
  #bin_boundaries=np.linspace(-100,100,n_bins+1)

  # bin width 0.05K
  n_samples=20000000 # 20,000,000
  #n_samples=50000000
  n_bins=20000
  bin_boundaries=np.linspace(-100,100,n_bins+1)

  # quick/dirty calcs
  if quick_and_dirty:
    n_samples=100000 # 100,000L
    n_bins=10000
    bin_boundaries=np.linspace(-1000,1000,n_bins+1)

  # priors sample for nice looking prior plots
  if prior_sample:
    # bin width 0.1K
    n_samples=1000000  #  1,000,000
    n_bins=2000
    bin_boundaries=np.linspace(-100,100,n_bins+1)

  # medium sample for more tractable calculations with toploop
  if medium_sample:
    # bin width 0.01K
    n_samples=1000000  #  1,000,000
    n_bins=20000
    bin_boundaries=np.linspace(-100,100,n_bins+1)

  # upped medium sample for more tractable calculations with toploop
  if upper_medium_sample:
    # bin width 0.01K
    n_samples=1500000  #  1,500,000
    n_bins=20000
    bin_boundaries=np.linspace(-100,100,n_bins+1)

  # small sample for more tractable calculations with toploop
  if small_sample:
    # bin width 0.05K
    n_samples=100000  #  100,000
    n_bins=20000
    bin_boundaries=np.linspace(-100,100,n_bins+1)

  # small sample for more tractable calculations with toploop
  if tiny_sample:
    # bin width 0.05K
    n_samples=300  #  300
    n_bins=20000
    bin_boundaries=np.linspace(-100,100,n_bins+1)

  #print ("bin_boundaries=",bin_boundaries)

  bin_width=bin_boundaries[1]-bin_boundaries[0]
  bin_centres=bin_boundaries[0:n_bins]+bin_width/2
  return (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)

def dnorm(x,mu,sd):
  exponent = -0.5*((x-mu)/sd)**2
  denominator = sd * (2*3.14159265359)**0.5
  return np.exp(exponent)/denominator

def dfat(x,mu,sd):
  # evalute student t PDF with 5 dof at x
  # note that as used hear mu is a vector of means and x is a single value
  # (following James Annan's R code)
  from scipy.stats import t
  print ("dfat type(x)=",type(x),x)
  print ("dfat type(mu)=",type(mu),mu)
  print ("dfat type(sd)=",type(sd),sd)
  return t.pdf((x-mu)/sd, 5)

def rnorm(nsamp,mean,sd):
  return np.random.normal(mean,sd,nsamp)

def rfat(nsamp,mean,sd):
  # return sample from student t distribution with 5 DoF
  sample=np.random.standard_t(5,nsamp)*sd+mean
  print ("rfat 5 DoF t-distribution mean and 5-95% range = ","{:10.2f}".format(sample.mean()),"{:10.2f}".format(sample.std()),"{:10.2f}".format(np.percentile(sample,5)),"{:10.2f}".format(np.percentile(sample,95)))
  return sample

def rt(nsamp,mean,sd,sample_size):
  sample=np.random.standard_t(sample_size-1,nsamp)*sd+mean

  print ("t-distribution mean and 5-95% range = ","{:10.2f}".format(sample.mean()),"{:10.2f}".format(sample.std()),"{:10.2f}".format(np.percentile(sample,5)),"{:10.2f}".format(np.percentile(sample,95)))

  return sample

def f2co2_process(normal01_2co2):

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  # values from Steve Klein:
  # updated 23/1/20
  fomu,fosigma=4.0,0.3
  f2co2_proc=fomu+normal01_2co2*fosigma
  #f2co2_proc=np.random.normal(fomu,fosigma,n_samples)

  return f2co2_proc

def normalise_pdf(pdf,bin_width):
  return(pdf/np.sum(pdf, dtype=np.float64)/bin_width)

def calc_uniform_s_prior(lower,upper,n_samples,weights):

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  s=np.random.uniform(lower,upper,n_samples)

  prior=sample2histogram(s,bin_boundaries,weights)

  return(prior)

#function to calculate the historical likelihood

def historical_dt_obs():

  DT_o = 1.03        # updated 19/3/20
  DT_sd = 0.14/1.64

  if alternative_hist_dt:
    DT_o=alternative_DT_o
    DT_sd=alternative_DT_sd

  print ("DT_o=",DT_o,"DT_sd=",DT_sd)

  return (DT_o,DT_sd)

def historical_likelihood(l_in,F2x,CO2_hist,F_other,historical_imbalance,lam_p):

  if calc_tecs:
    # use transfer function if calculating true equilibrium climate sensitivity (TECS)
    # in this case prior l_in is lambda_TECS
    # TECS = S * (1+transfer)
    # lambda_TECS = lambda_S / (1+transfer)
    # lambda_S = lambda_TECS * (1+transfer)
    # historical evidence is a function of l which is lambda_S
    # so we want l = lambda_S = l_in*(1+transfer)
    l=l_in*(1+transfer)
  else:
    l=l_in

  (DT_o,DT_sd) = historical_dt_obs()

  if reduce_dt_dn_uncert:
    #DT_sd=DT_sd/2.0
    #mean_historical_imbalance=np.mean(historical_imbalance)
    #historical_imbalance=mean_historical_imbalance+(historical_imbalance-mean_historical_imbalance)/2
    DT_sd=DT_sd/10.0
    mean_historical_imbalance=np.mean(historical_imbalance)
    historical_imbalance=mean_historical_imbalance+(historical_imbalance-mean_historical_imbalance)/10

  #model for temperature change
  DT = -1*(CO2_hist*F2x/np.mean(F2x) + F_other - historical_imbalance)/(l + lam_p)

  if reduce_fhist_uncert:
    hist_forc=CO2_hist*F2x/np.mean(F2x) + F_other
    mean_hist_forc=np.mean(hist_forc)
    hist_forc=mean_hist_forc+(hist_forc-mean_hist_forc)/10
    DT = -1*(hist_forc - historical_imbalance)/(l + lam_p)

  if no_fhist_uncert:
    hist_forc=CO2_hist*F2x/np.mean(F2x) + F_other
    hist_forc=np.mean(hist_forc)
    DT = -1*(hist_forc - historical_imbalance)/(l + lam_p)

  #likelihood based on observed temperature rise of 0.91 with uncertainty of 0.14 at 90% = 1.64sd of gaussian

  likelihood=dnorm(DT_o,DT,DT_sd)

  return(likelihood)

def historical_classic_sample(F2x,CO2_hist,F_other,historical_imbalance,lam_p):

  (DT_o,DT_sd) = historical_dt_obs()

  #model for temperature change
  DT=np.random.normal(DT_o,DT_sd,n_samples)

  # DT = -1*(CO2_hist*F2x/np.mean(F2x) + F_other - historical_imbalance)/(l + lam_p)
  # (l + lam_p) =  -1*(CO2_hist*F2x/np.mean(F2x) + F_other - historical_imbalance)/DT
  l =  -1*(CO2_hist*F2x/np.mean(F2x) + F_other - historical_imbalance)/DT - lam_p
  s = -1 * F2x / l

  return(s,l)

if plot_priors:

  # calculate priors only

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()
  equal_weights=np.full(n_samples,1.0)
  normal01_2co2=np.random.normal(0,1,n_samples)
  uprior=calc_uniform_s_prior(0,20,n_samples,equal_weights)

  #90% (5-95%) confidence interval for normal distribution is +/-1.645 * sigma

  #calculate gaussian samples to correlate CO2 forcings

  ticklen=0.05
  xrange=[-1,21]
  yrange=[0,1]

  ul_prior_label='UL Prior U(-20,0)'
  (prior,f2co2_proc,sampled_l_prior,sampled_s_prior)=calc_prior_uniform_lambda(n_samples,normal01_2co2,equal_weights)

  f_prior_pdf , bin_boundaries = np.histogram(f2co2_proc, bins=bin_boundaries)
  #f_prior_pdf = normalise_pdf(f_prior_pdf,bin_width)

  s_prior_pdf_ul , bin_boundaries = np.histogram(sampled_s_prior, bins=bin_boundaries)
  #s_prior_pdf_ul = normalise_pdf(s_prior_pdf_ul,bin_width)

  l_prior_pdf_ul , bin_boundaries = np.histogram(sampled_l_prior, bins=bin_boundaries)
  #l_prior_pdf_ul = normalise_pdf(l_prior_pdf_ul,bin_width)

  uli_prior_label='ULI Priors ~ U(-10,-10)'
  (prior,f2co2_proc,sampled_l_prior,sampled_s_prior,l1,l2,l3,l4,l5,l6)=calc_prior_uniform_lambda_is(n_samples,normal01_2co2,equal_weights,False)

  s_prior_pdf_uli , bin_boundaries = np.histogram(sampled_s_prior, bins=bin_boundaries)
  #s_prior_pdf_uli = normalise_pdf(s_prior_pdf_uli,bin_width)

  l_prior_pdf_uli , bin_boundaries = np.histogram(sampled_l_prior, bins=bin_boundaries)
  #l_prior_pdf_uli = normalise_pdf(l_prior_pdf_uli,bin_width)

  process_bu_prior_label='Process BU Prior L~N(-1.35,0.41)'
  (prior,f2co2_proc,sampled_l_prior,sampled_s_prior)=calc_process_bu_prior(n_samples,normal01_2co2,equal_weights)

  s_prior_pdf_process_bu , bin_boundaries = np.histogram(sampled_s_prior, bins=bin_boundaries)
  #s_prior_pdf_process_bu = normalise_pdf(s_prior_pdf_process_bu,bin_width)

  l_prior_pdf_process_bu , bin_boundaries = np.histogram(sampled_l_prior, bins=bin_boundaries)
  #l_prior_pdf_process_bu = normalise_pdf(l_prior_pdf_process_bu,bin_width)

  us_uli_prior_label='US2 Prior (0,20)'
  (prior,f2co2_proc,l1,l2,l3,l4,l5,l6,sampled_s_prior,prior_weights) = \
    calc_uniform_s_prior_from_li(n_samples,normal01_2co2,equal_weights,False,'ULI',weight_prior)
  sampled_l_prior=l1+l2+l3+l4+l5+l6

  s_prior_pdf_us_uli , bin_boundaries = np.histogram(sampled_s_prior, bins=bin_boundaries,weights=prior_weights)
  #s_prior_pdf_us_uli = normalise_pdf(s_prior_pdf_us_uli,bin_width)

  l_prior_pdf_us_uli , bin_boundaries = np.histogram(sampled_l_prior, bins=bin_boundaries,weights=prior_weights)
  #l_prior_pdf_us_uli = normalise_pdf(l_prior_pdf_us_uli,bin_width)

  f_prior_pdf_us_uli , bin_boundaries = np.histogram(f2co2_proc, bins=bin_boundaries)
  #f_prior_pdf_us_uli = normalise_pdf(f_prior_pdf_us_uli,bin_width)

  us_prior_label='US1 Prior U(0,20)'
  (prior,f2co2_proc,sampled_l_prior,sampled_s_prior)=calc_prior_uniform_s(n_samples,normal01_2co2,equal_weights)

  s_prior_pdf_us , bin_boundaries = np.histogram(sampled_s_prior, bins=bin_boundaries)
  #s_prior_pdf_us = normalise_pdf(s_prior_pdf_us,bin_width)

  l_prior_pdf_us , bin_boundaries = np.histogram(sampled_l_prior, bins=bin_boundaries)
  #l_prior_pdf_us = normalise_pdf(l_prior_pdf_us,bin_width)

if plot_priors:

  # calculate and plot priors only

  plt.close('all')

  legend_spacing=0.08
  # Two subplots, the axes array is 1-d
  #f, axarr = plt.subplots(1,2, sharey=True)

  color='k'
  #axarr[0]#plt.plot(
  #plt.plot([4.5,4.5],yrange,color=color)

  #axarr[0].set_xlim([-6,2])
  #axarr[0].set_ylim([0,1])
  #axarr[0].set_xlabel('L (Wm-2K-1)')

  yrange=[0,2]
  plt.xlim([-6,2])
  plt.ylim(yrange)
  plt.xlabel('L (Wm-2K-1)')

  # **** plot l prior predictive distributions *****

  xpos=-5.5
  pos=2
  dy=1
  label='UL'
  color='r'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),ul_prior_label,color=color)
  plt.plot(bin_centres,l_prior_pdf_ul,color)

  pos=pos+dy
  label='ULI'
  color='b'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),uli_prior_label,color=color)
  plt.plot(bin_centres,l_prior_pdf_uli,color)

  pos=pos+dy
  label='Process BU'
  color='tab:orange'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),process_bu_prior_label,color=color)
  plt.plot(bin_centres,l_prior_pdf_process_bu,color)

  pos=pos+dy
  label='US1'
  color='c'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),us_prior_label,color=color)
  plt.plot(bin_centres,l_prior_pdf_us,color,linewidth=5)

  #plt.plot([1.5,1.5],yrange,color='k')
  #plt.plot([4.5,4.5],yrange,color='k')

  pos=pos+dy
  label='US2'
  color='k'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),us_uli_prior_label,color=color)
  plt.plot(bin_centres,l_prior_pdf_us_uli,color,linewidth=3)

  xpos=pos+dy

  plt.savefig(outpath+'prior_plot_l.png')
  plt.close('all')

  # **** plot S prior predictive distributions *****

  plt.ylim(yrange)
  plt.xlim([-2,21])
  plt.xlabel('S (K)')
  xpos=2
  pos=2
  label='UL'
  color='r'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),ul_prior_label,color=color)
  plt.plot(bin_centres,s_prior_pdf_ul,color)

  pos=pos+dy
  label='ULI'
  color='b'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),uli_prior_label,color=color)
  plt.plot(bin_centres,s_prior_pdf_uli,color)

  pos=pos+dy
  label='Process BU'
  color='tab:orange'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),process_bu_prior_label,color=color)
  plt.plot(bin_centres,s_prior_pdf_process_bu,color)

  pos=pos+dy
  label='US1'
  color='c'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),us_prior_label,color=color)
  plt.plot(bin_centres,s_prior_pdf_us,color,linewidth=5)

  pos=pos+dy
  label='US2'
  color='k'
  plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),us_uli_prior_label,color=color)
  plt.plot(bin_centres,s_prior_pdf_us_uli,color,linewidth=3)

  #pos=pos+dy
  #label='F_US_ULI (-)'
  #color='c'
  #plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),label,color=color)
  #plt.plot(bin_centres,f_prior_pdf_us_uli,color,linestyle='-')

  #pos=pos+dy
  #label='F2CO2 (-)'
  #color='r'
  #plt.text(xpos,yrange[1]*1-legend_spacing*(pos+1),label,color=color)
  #plt.plot(bin_centres,f_prior_pdf,color,linestyle='-')

  print('XXprior')

  plt.savefig(outpath+'prior_plot_s.png')
  plt.close('all')

def likelihood(uniform_sample,weights,label,var):

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  likelihood , bin_boundaries = np.histogram(uniform_sample, weights=weights, bins=bin_boundaries)

  # normalise to max in the plotting range
  #tmp=np.copy(likelihood)
  #tmp[bin_centres > 10]=0
  #tmp[bin_centres < -10]=0
  #likelihood=likelihood/tmp.max()

  stats=1
  if stats:
    #print ("likelihood.shape=",likelihood.shape)
    #print ("bin_centres.shape=",bin_centres.shape)
    #print ("bin_boundaries.shape=",bin_boundaries.shape)

    peak_distance=(likelihood-1.0)**2
    peak=bin_centres[peak_distance == peak_distance.min()]

    peak=np.mean(peak)
    print ("peak=",peak)

    lower_half=likelihood[bin_centres < peak]
    lower_bins=bin_centres[bin_centres < peak]
    distance=(lower_half-0.2)**2
    lower_crossing=lower_bins[distance == distance.min()]
    lower_crossing=lower_crossing[0]

    upper_half=likelihood[bin_centres > peak]
    upper_bins=bin_centres[bin_centres > peak]
    distance=(upper_half-0.2)**2
    upper_crossing=upper_bins[distance == distance.min()]
    upper_crossing=upper_crossing[0]

    print ("True Likelihood for",var,label,lower_crossing,peak,upper_crossing)

  return (likelihood)

#from scipy.interpolate import UnivariateSpline

#from statsmodels.nonparametric.smoothers_lowess import lowess

#from scipy.signal import savgol_filter

def pseudo_likelihood(updated_pdf,prior_pdf,label,var,smooth=True):

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

  pseudo_likelihood = updated_pdf / prior_pdf

  window=41
  order=2
  #s1 = savgol_filter(updated_pdf, window, 2)
  #s2 = savgol_filter(prior_pdf, window, 2)
  #pseudo_likelihood = s1/s2

  #k=5
  #s1=UnivariateSpline(bin_centres, updated_pdf, s=k)
  #s2=UnivariateSpline(bin_centres, prior_pdf, s=k)

  #pseudo_likelihood=s1(bin_centres)/s2(bin_centres)

  pseudo_likelihood[np.isnan(pseudo_likelihood)] = 0.

  #k=5
  #s=UnivariateSpline(bin_centres, pseudo_likelihood, s=k)
  #pseudo_likelihood=s(bin_centres)

  #pseudo_likelihood = lowess(pseudo_likelihood, bin_centres, is_sorted=True, frac=0.025, it=0)
  #pseudo_likelihood = lowess(pseudo_likelihood, bin_centres)

  #if smooth:
  #  pseudo_likelihood = savgol_filter(pseudo_likelihood, window, order)

  # normalise to max in the plotting range
  #tmp=np.copy(pseudo_likelihood)
  #tmp[bin_centres > 10]=0
  #tmp[bin_centres < -10]=0
  #pseudo_likelihood=pseudo_likelihood/tmp.max()

  stats=0
  if stats:
    peak_distance=(pseudo_likelihood-1.0)**2
    peak=bin_centres[peak_distance == peak_distance.min()]

    peak=np.mean(peak)
    print ("peak=",peak)

    lower_half=pseudo_likelihood[bin_centres < peak]
    lower_bins=bin_centres[bin_centres < peak]
    distance=(lower_half-0.2)**2
    lower_crossing=lower_bins[distance == distance.min()]
    lower_crossing=lower_crossing[0]

    upper_half=pseudo_likelihood[bin_centres > peak]
    upper_bins=bin_centres[bin_centres > peak]
    distance=(upper_half-0.2)**2
    upper_crossing=upper_bins[distance == distance.min()]
    upper_crossing=upper_crossing[0]

    print ("Pseudo likelihood for",var,label,lower_crossing,peak,upper_crossing)

  return (pseudo_likelihood)


def plot_prior_likelihoods_posterior():

  # plot prior, likelihoods and posterior every

  import matplotlib.pyplot as plt

  plt.rcParams['figure.figsize'] = [10, 7]

  plt.close('all')

  biplot=0

  if biplot:
    # Two subplots, the axes array is 1-d
    f, axarr = plt.subplots(1,2, sharey=True)

  color='k'

  # **** lambda plot ******

  l_prior=np.copy(full_l_prior_pdf)
  l_prior_label='Lambda '+prior_label

  xpos=-11
  pos=0.8
  dy=0.07


  color='k'
  #plt.plot([1.5,1.5],yrange,color=color)
  #plt.plot([4.5,4.5],yrange,color=color)

  #pos=2
  #color='k'
  #label='s_prior_pdf'
  #plt.text(-1,1-0.05*(pos+1),label,color='r')
  color='r'
  #plt.plot(bin_centres,s_prior_pdf,color)

  if biplot:
    plt=axarr[0]
    plt.set_xlim([-6,1])
    plt.set_ylim([0,1])
    plt.set_xlabel('L (Wm-2K-1)')
  else:
    plt.xlim([-6,1])
    plt.ylim([0,1])
    plt.xlabel('L (Wm-2K-1)')

  if plot_prior:
    plt.text(xpos,yrange[1]*1-0.05*(pos+1),l_prior_label,color=color)
    plt.plot(bin_centres,l_prior_pdf,color)
    plt.plot(bin_centres,full_l_prior_pdf,color,linestyle='dotted')

  pos=3

  if not no_process_bu:
    pos=pos+1
    color='g'
    label='BU Process Likelihood'
    plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
    plt.plot(bin_centres,l_process_bu_likelihood,color,linewidth=2)
  else:
    pos=pos+0.5

  if not no_process_ec:
    pos=pos+1
    color='y'
    label='EC Process Likelihood'
    plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
    plt.plot(bin_centres,l_process_ec_likelihood,color,linewidth=2)
  else:
    pos=pos+0.5

  pos=pos+1
  color='b'
  label='Historical Likelihood'
  plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
  plt.plot(bin_centres,l_hist_likelihood,color,linewidth=2)

  pos=pos+1
  color='c'
  label='Cold Paleo Likelihood'
  plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
  plt.plot(bin_centres,l_paleo_cold_likelihood,color,linewidth=2)

  pos=pos+1
  color='m'
  label='Hot Paleo Likelihood'
  plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
  plt.plot(bin_centres,l_paleo_hot_likelihood,color,linewidth=2)

  color='r'
  print ('statistics for l_prior:')

  if plot_prior:
    #plt.plot(bin_centres,l_prior,linewidth=2,color=color)
    #plt.text(xpos,yrange[1]*0.9,l_prior_label,color=color)

    ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(full_l_prior_pdf,bin_boundaries,bin_width,'full_l_prior_pdf')

    #plt.text(xpos,yrange[1]*0.85,'P(ECS < 1.5) = '+"{:.2f}".format(p_ecs_lt_1p5/100.0),color=color)
    #plt.text(xpos,yrange[1]*0.8,'P(ECS > 4.5) = '+"{:.2f}".format(p_ecs_gt_4p5/100.0),color=color)
    plt.text(xpos,yrange[1]*0.85,"L 5th %ile = "+"{:.1f}".format(ecs5),color=color)
    plt.text(xpos,yrange[1]*0.8,'L 95th %ile = '+"{:.1f}".format(ecs95),color=color)

  print ( "statistics for l_posterior:" )

  if plot_posterior:
    color='k'
    plt.plot(bin_centres,l_posterior,linewidth=2,color=color)
    plt.text(xpos,yrange[1]*0.5,l_posterior_label,color=color)

    ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(l_posterior,bin_boundaries,bin_width,'L Posterior')

    #plt.text(xpos,yrange[1]*0.45,'L 5th %ile = '+"{:.1f}".format(ecs5),color=color)
    #plt.text(xpos,yrange[1]*0.4,'L 95th %ile = '+"{:.1f}".format(ecs95),color=color)

    print ('xx l posterior',toploop_index)
    #if biplot == 0:
  plt.savefig(out_dir+'/'+calc_id+'2.png')
  plt.close('all')

  # **** S plot ******

  if biplot:
    plt=axarr[1]
    plt.set_xlim([-1,21])
    plt.set_ylim([0,1])
    plt.set_xlabel('S (K)')
  else:
    plt.xlim([-1,21])
    plt.ylim([0,1])
    plt.xlabel('S (K)')

  xpos=8
  pos=3

  if not no_process_bu:
    pos=pos+1
    color='g'
    label='BU Process Likelihood'
    plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
    plt.plot(bin_centres,s_process_bu_likelihood,color,linewidth=2)
  else:
    pos=pos+0.5

  if not no_process_ec:
    pos=pos+1
    color='y'
    label='EC Process Likelihood'
    plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
    plt.plot(bin_centres,s_process_ec_likelihood,color,linewidth=2)
  else:
    pos=pos+0.5

  pos=pos+1
  color='b'
  label='Historical Likelihood'
  plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
  plt.plot(bin_centres,s_hist_likelihood,color,linewidth=2)

  pos=pos+1
  color='c'
  label='Cold Paleo Likelihood'
  plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
  plt.plot(bin_centres,s_paleo_cold_likelihood,color,linewidth=2)

  pos=pos+1
  color='m'
  label='Hot Paleo Likelihood'
  plt.text(xpos,yrange[1]*1-0.05*(pos+1),label,color=color)
  plt.plot(bin_centres,s_paleo_hot_likelihood,color,linewidth=2)

  color='r'
  print ('statistics for prior:')

  if plot_prior:
    plt.plot(bin_centres,s_prior_pdf,linewidth=2,color=color)
    plt.plot(bin_centres,full_s_prior_pdf,color,linestyle='dotted')

    plt.text(xpos,yrange[1]*0.9,prior_label,color=color)

    ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(s_prior_pdf,bin_boundaries,bin_width,'s_prior_pdf')

    ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(full_s_prior_pdf,bin_boundaries,bin_width,'full_s_prior_pdf')

    #plt.text(xpos,yrange[1]*0.85,'P(ECS < 1.5) = '+"{:.2f}".format(p_ecs_lt_1p5/100.0),color=color)
    #plt.text(xpos,yrange[1]*0.8,'P(ECS > 4.5) = '+"{:.2f}".format(p_ecs_gt_4p5/100.0),color=color)
    plt.text(xpos,yrange[1]*0.85,'ECS 5th %ile = '+"{:.1f}".format(ecs5),color=color)
    plt.text(xpos,yrange[1]*0.8,'ECS 95th %ile = '+"{:.1f}".format(ecs95),color=color)

  print ('statistics for posterior:')

  if plot_posterior:
    color='k'
    plt.plot(bin_centres,posterior,linewidth=2,color=color)
    plt.text(xpos,yrange[1]*0.5,s_posterior_label,color=color)

    ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(posterior,bin_boundaries,bin_width,'Not Sure Posterior')

    #plt.text(xpos,yrange[1]*0.45,'P(ECS < 1.5) = '+"{:.2f}".format(p_ecs_lt_1p5/100.0),color=color)
    #plt.text(xpos,yrange[1]*0.4,'P(ECS > 4.5) = '+"{:.2f}".format(p_ecs_gt_4p5/100.0),color=color)
    plt.text(xpos,yrange[1]*0.45,'ECS 5th %ile = '+"{:.1f}".format(ecs5),color=color)
    plt.text(xpos,yrange[1]*0.4,'ECS 95th %ile = '+"{:.1f}".format(ecs95),color=color)

  print ('xx s posterior',toploop_index)
  #if biplot == 0:
  plt.savefig(out_dir+'/'+calc_id+'3.png')
  plt.close('all')

if plot_ecs_pdf:

  (n_bins,bin_boundaries,bin_centres,bin_width,n_samples) = setup_bins()

  (df_recent,df_baseline,uncs_neg,uncs_pos) = read_hist_forc()

  full_importance_density_scaling = np.array([0])

  for toploop_index in range(0,ntop):

    # ----------------------------------------------------------------------------
    # set up some widely used variables
    # ----------------------------------------------------------------------------

    print ('xx toploop_index=', toploop_index)
    print (datetime.datetime.now())

    (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()
    equal_weights=np.full(n_samples,1.0)
    normal01_2co2=np.random.normal(0,1,n_samples)
    norm2=np.random.normal(0,1,n_samples)

    # ----------------------------------------------------------------------------
    # calculate prior
    # ----------------------------------------------------------------------------

    #uprior=calc_uniform_s_prior(0,20,n_samples)

    #90% (5-95%) confidence interval for normal distribution is +/-1.645 * sigma

    # calculate gaussian samples to correlate CO2 forcings

    ticklen=0.05
    xrange=[-1,21]
    yrange=[0,1]
    full_prior_weights=np.copy(equal_weights)

    if prior_id == 'Process':
      prior_label='Process BU Prior'
      (prior,f2co2_proc,sampled_l_prior,sampled_s_prior)=calc_process_prior(n_samples,normal01_2co2,equal_weights)
      full_sampled_l_prior=np.copy(sampled_l_prior)
      full_sampled_s_prior=np.copy(sampled_s_prior)
      prior_weights=np.copy(equal_weights)

    if prior_id == 'UL':
      prior_label='UL Prior [-20,0]'
      (prior,f2co2_proc,sampled_l_prior,sampled_s_prior)=calc_prior_uniform_lambda(n_samples,normal01_2co2,equal_weights)
      full_sampled_l_prior=np.copy(sampled_l_prior)
      full_sampled_s_prior=np.copy(sampled_s_prior)
      prior_weights=np.copy(equal_weights)

    if prior_id == 'PROCBU':
      prior_label='L Prior ~N(-1.35,0.41)'
      (prior,f2co2_proc,sampled_l_prior,sampled_s_prior)=calc_process_bu_prior(n_samples,normal01_2co2,equal_weights)
      full_sampled_l_prior=np.copy(sampled_l_prior)
      full_sampled_s_prior=np.copy(sampled_s_prior)
      prior_weights=np.copy(equal_weights)

    if prior_id == 'US1':
      prior_label='US1 Prior [0,20]'
      (prior,f2co2_proc,sampled_l_prior,sampled_s_prior)=calc_prior_uniform_s(n_samples,normal01_2co2,equal_weights)
      full_sampled_l_prior=np.copy(sampled_l_prior)
      full_sampled_s_prior=np.copy(sampled_s_prior)
      prior_weights=np.copy(equal_weights)

    if prior_id == 'ULI':
      prior_label='ULI Priors [-10,10]'
      (prior,f2co2_proc,sampled_l_prior,sampled_s_prior,l1,l2,l3,l4,l5,l6)=\
        calc_prior_uniform_lambda_is(n_samples,normal01_2co2,equal_weights,restrict_li_range)
      if restrict_li_range:
        # Calculate full prior to show also
        (full_s_prior_pdf,f2co2_proc,full_sampled_l_prior,full_sampled_s_prior,
         full_l1,full_l2,full_l3,full_l4,full_l5,full_l6) = \
         calc_prior_uniform_lambda_is(n_samples,normal01_2co2,equal_weights,False)
        print ('prior_id=',prior_id)
        print ('calc_id=',calc_id)
        print ('full_sampled_s_prior.shape[0]=',full_sampled_s_prior.shape[0])
      else:
        full_sampled_l_prior=np.copy(sampled_l_prior)
        full_sampled_s_prior=np.copy(sampled_s_prior)
      prior_weights = np.copy(equal_weights)

    if prior_id == 'US2':
      prior_label='US2 Prior [0,20]'
      # restrict prior to process evidence neighbourhood to make 6d calculation tractable
      if restrict_li_range:
        # Calculate full prior to pass into restricted li calculation
        (full_prior,f2co2_proc,full_l1,full_l2,full_l3,full_l4,full_l5,full_l6,full_sampled_s_prior,full_prior_weights,full_importance_density_scaling) = \
        calc_uniform_s_prior_from_li(n_samples,normal01_2co2,equal_weights,False,input_prior,weight_prior,full_importance_density_scaling)
        full_sampled_l_prior=full_l1+full_l2+full_l3+full_l4+full_l5+full_l6

        # Do restricted li calculation
        (prior,f2co2_proc,l1,l2,l3,l4,l5,l6,sampled_s_prior,prior_weights,full_importance_density_scaling) = \
          calc_uniform_s_prior_from_li(n_samples,normal01_2co2,equal_weights,restrict_li_range,input_prior,weight_prior,full_importance_density_scaling)
        sampled_l_prior=l1+l2+l3+l4+l5+l6
      else:
        (prior,f2co2_proc,l1,l2,l3,l4,l5,l6,sampled_s_prior,prior_weights,full_importance_density_scaling) = \
          calc_uniform_s_prior_from_li(n_samples,normal01_2co2,equal_weights,restrict_li_range,input_prior,weight_prior,full_importance_density_scaling)
        sampled_l_prior=l1+l2+l3+l4+l5+l6

        full_prior=np.copy(prior)
        full_sampled_l_prior=np.copy(sampled_l_prior)
        full_prior_weights=np.copy(prior_weights)
        full_l1=np.copy(l1)
        full_l2=np.copy(l2)
        full_l3=np.copy(l3)
        full_l4=np.copy(l4)
        full_l5=np.copy(l5)
        full_l6=np.copy(l6)
        full_sampled_s_prior=np.copy(sampled_s_prior)


      print ("correlation l1 l2=",np.corrcoef(full_l1,full_l2))
      print ("correlation l2 l3=",np.corrcoef(full_l2,full_l3))
      print ("correlation l3 l4=",np.corrcoef(full_l3,full_l4))
      print ("correlation l4 l5=",np.corrcoef(full_l4,full_l5))
      print ("correlation l1 l2+l3+l4+l5+l6=",np.corrcoef(full_l1,full_l2+full_l3+full_l4+full_l5+full_l6))
      s=np.copy(full_sampled_s_prior)

      low=2.9
      high=3.1
      sub1=full_l1[s<high]
      sub2=full_l2[s<high]
      sub3=full_l3[s<high]
      sub4=full_l4[s<high]
      sub5=full_l5[s<high]
      sub6=full_l6[s<high]
      subs=s[s<high]
      sub1=np.array(sub1[subs>low])
      sub2=np.array(sub2[subs>low])
      sub3=np.array(sub3[subs>low])
      sub4=np.array(sub4[subs>low])
      sub5=np.array(sub5[subs>low])
      sub6=np.array(sub6[subs>low])

      #print ("sub5=",sub5)
      #print ("sub5.shape()=",sub5.shape())

      print ("correlation l1 l2|S~=",np.corrcoef(sub1,sub2))
      print ("correlation l2 l3|S~=3:",np.corrcoef(sub2,sub3))
      print ("correlation l3 l4|S~=3:",np.corrcoef(sub3,sub4))
      print ("correlation l4 l5|S~=3:",np.corrcoef(sub4,sub5))
      print ("correlation l1 l2+l3+l4+l5+l6|S~=3:",np.corrcoef(sub1,sub2+sub3+sub4+sub5+sub6))

      low=8.9
      high=9.1
      sub1=full_l1[s<high]
      sub2=full_l2[s<high]
      sub3=full_l3[s<high]
      sub4=full_l4[s<high]
      sub5=full_l5[s<high]
      sub6=full_l6[s<high]
      subs=s[s<high]
      sub1=np.array(sub1[subs>low])
      sub2=np.array(sub2[subs>low])
      sub3=np.array(sub3[subs>low])
      sub4=np.array(sub4[subs>low])
      sub5=np.array(sub5[subs>low])
      sub6=np.array(sub6[subs>low])

      #print ("sub5=",sub5)
      #print ("sub5.shape()=",sub5.shape())

      print ("correlation l1 l2|S~=9",np.corrcoef(sub1,sub2))
      print ("correlation l2 l3|S~=9:",np.corrcoef(sub2,sub3))
      print ("correlation l3 l4|S~=9:",np.corrcoef(sub3,sub4))
      print ("correlation l4 l5|S~=9:",np.corrcoef(sub4,sub5))
      print ("correlation l1 l2+l3+l4+l5+l6|S~=9:",np.corrcoef(sub1,sub2+sub3+sub4+sub5+sub6))

    print ('prior_id=',prior_id)
    print ('calc_id=',calc_id)
    print ('full_sampled_s_prior.shape[0]=',full_sampled_s_prior.shape[0])
    print ('full_prior_weights=',full_prior_weights.shape[0])
    print ('sampled_s_prior.shape[0]=',sampled_s_prior.shape[0])
    assert (full_sampled_s_prior.shape[0] == full_prior_weights.shape[0])
    full_s_prior_pdf=sample2histogram(full_sampled_s_prior,bin_boundaries,full_prior_weights)
    assert (sampled_s_prior.shape[0] == full_prior_weights.shape[0])
    s_prior_pdf = sample2histogram(sampled_s_prior,bin_boundaries,prior_weights)

    s_pdf=np.copy(s_prior_pdf)
    s_sample=np.copy(sampled_s_prior)

    l_prior_pdf = sample2histogram(sampled_l_prior, bin_boundaries,prior_weights)
    full_l_prior_pdf = sample2histogram(full_sampled_l_prior, bin_boundaries,full_prior_weights)

    l_prior = np.copy(full_l_prior_pdf)

    l_pdf=np.copy(l_prior_pdf)
    l_sample=np.copy(sampled_l_prior)

    # ----------------------------------------------------------------------------
    # Set up transfer function between S and true ecs
    # ----------------------------------------------------------------------------
    #this is the uncertain transfer function between the true equilibrium climate sensitivity
    #and the regression value S that we target
    #implemented as true ECS = regression * (1+transfer) based on Maria Rugenstein's work
    #hence ECS / S ~ N(1.06,0.2)

    # If we are calculating S this us used with the paleo evidence to convert
    # likelihoods for equilibrium paleo sensitivity to likelihoods for S
    # If we are calculating true ECS this is used with the historical and process evidence
    # to convert likelihoods for S to likelihoods for true ECS

    nsamp = n_samples
    if fat_tails:
      transfer = rfat(nsamp,mean=0.06,sd=0.2)
    else:
      transfer = rnorm(nsamp,mean=0.06,sd=0.2)

    # calculate weighted and unweighted prior PDFs for the transfer function
    # use f bins because need precision over a small range

    f_bin_boundaries=np.linspace(-10,10,n_bins+1)
    f_bin_width=f_bin_boundaries[1]-f_bin_boundaries[0]
    f_bin_centres=f_bin_boundaries[0:n_bins]+f_bin_width/2

    transfer_weighted_prior_pdf, f_bin_boundaries = np.histogram(transfer, weights=full_prior_weights, bins=f_bin_boundaries)
    transfer_unweighted_prior_pdf, f_bin_boundaries = np.histogram(transfer, bins=f_bin_boundaries)

    # ----------------------------------------------------------------------------
    # Prior sample for ECS
    # ----------------------------------------------------------------------------

    ecs_sample=s_sample * (1+transfer)

    # ----------------------------------------------------------------------------
    # Set up variables needed for likelihood calculations
    # ----------------------------------------------------------------------------

    (n_bins,bin_boundaries,bin_centres,bin_width,n_samples)=setup_bins()

    if weight_prior:
      weights=np.copy(prior_weights)
    else:
      weights=np.copy(equal_weights)

    li_prior=False
    if prior_id == 'ULI' or prior_id == 'US2':
      li_prior=True

    # skip calculating marginal process likelihoods for lambda components.
    # We can't do this until we have code to restrict the ranges of the li
    # components for all priors.
    # The alternatice is to se pseudo-likelihoods
    calculate_li_marginal_process_likelihoods=0

    # set up uniform samples for calculation of likelihoods
    if calculate_li_marginal_process_likelihoods and li_prior:
      # Sample f2co2 and lambda_i components which imply a uniform distribution in l for use in calculating likelihoods
      (tmp_prior,tmp_f2co2_proc,ul_l1,ul_l2,ul_l3,ul_l4,ul_l5,ul_l6,sampled_s_uniform_l) = calc_uniform_l_prior_from_li(n_samples,normal01_2co2,equal_weights)
      sampled_l_uniform_l=ul_l1+ul_l2+ul_l3+ul_l4+ul_l5+ul_l6
      plot_sample(sampled_l_uniform_l,'should be uniform in l')
      # Sample f2co2 and lambda_i components which imply a uniform distribution in S *not a prior*
      (tmp_prior,tmp_f2co2_proc,us_l1,us_l2,us_l3,us_l4,us_l5,ul_l6,sampled_s_uniform_s) = calc_uniform_s_prior_from_li(n_samples,normal01_2co2,equal_weights)
      sampled_l_uniform_s=us_l1+us_l2+us_l3+us_l4+us_l5+us_l6
    else:
      (uniform_s_density_for_s,f2co2_proc,sampled_l_uniform_s,sampled_s_uniform_s)=calc_prior_uniform_s(n_samples,normal01_2co2,equal_weights)
      (uniform_l_density_for_s,f2co2_proc,sampled_l_uniform_l,sampled_s_uniform_l)=calc_prior_uniform_lambda(n_samples,normal01_2co2,equal_weights)

    # ----------------------------------------------------------------------------
    # process update
    # ----------------------------------------------------------------------------

    #f2co2_proc=f2co2_process(normal01_2co2)

    if no_process_bu:
      print ('No process BU likelihood with process prior')
      s_process_bu_likelihood=np.full(n_bins-1,0.0)
      l_process_bu_likelihood=np.full(n_bins-1,0.0)
    else:
      print ('Applying BU likelihood')
      # bayesian update for BU process evidence

      if li_prior:
        print ('xx bayesian update for BU process evidence')
        print ("XX2 l5.shape=",l5.shape)
        process_likelihood_bu_weights=process_likelihood_weight_li(l1,l2,l3,l4,l5,l6,plot_li_likelihoods)
      else:
        if pattern_effect_dependence2:
          process_likelihood_bu_weights=process_likelihood_weight_l_pattern_effect(sampled_l_prior,lam_p_p)
        else:
          process_likelihood_bu_weights=process_likelihood_weight_l(sampled_l_prior)

      weights = weights * process_likelihood_bu_weights

      s_pdf_updated , bin_boundaries = np.histogram(sampled_s_prior, weights=weights, bins=bin_boundaries)
      #s_pdf_updated = normalise_pdf(s_pdf_updated,bin_width)

      s_pdf=np.copy(s_pdf_updated)

      l_pdf_updated , bin_boundaries = np.histogram(sampled_l_prior, weights=weights, bins=bin_boundaries)
      #l_pdf_updated = normalise_pdf(l_pdf_updated,bin_width)

      print ('Calculate process pseudo likelihood for L')
      l_process_bu_likelihood = pseudo_likelihood(l_pdf_updated,l_pdf,'Process BU','L')

      l_pdf=np.copy(l_pdf_updated)

      print ('Calculate process pseudo likelihood for S')
      s_pdf_prior_plus_process_bu = sample2histogram(sampled_s_prior, bin_boundaries, weights=process_likelihood_bu_weights)
      s_process_bu_likelihood = pseudo_likelihood(s_pdf_prior_plus_process_bu,s_prior_pdf,'Process BU','S')

      #if li_prior:
      #  if calculate_li_marginal_process_likelihoods:
      #    process_likelihood_bu_weights=process_likelihood_weight_li(us_l1,us_l2,us_l3,us_l4,us_l5,us_l6,plot_li_likelihoods)
      #    s_process_bu_likelihood=likelihood(sampled_s_uniform_s,process_likelihood_bu_weights,'Process BU','S')
      #else:
      #  process_likelihood_bu_weights=process_likelihood_weight_l(sampled_l_uniform_s)
      #  s_process_bu_likelihood=likelihood(sampled_s_uniform_s,process_likelihood_bu_weights,'Process BU','S')

      # calculate pseudo likelihood using process_likelihood_bu_weights based on prior from bayesian update

      #print ('Calculate process likelihood for l')

      #if li_prior:
      #  if calculate_li_marginal_process_likelihoods:
      #    process_likelihood_bu_weights=process_likelihood_weight_li(ul_l1,ul_l2,ul_l3,ul_l4,ul_l5,ul_l6,plot_li_likelihoods)
      #    l_process_bu_likelihood=likelihood(sampled_l_uniform_l,process_likelihood_bu_weights,'Process BU','L')
      #  else:
      #    # calculate pseudo likelihood using process_likelihood_bu_weights based on prior from bayesian update
      #    l_pdf_prior_plus_process_bu = sample2histogram(sampled_l_prior, bin_boundaries, weights=process_likelihood_bu_weights)
      #    l_process_bu_likelihood = pseudo_likelihood(l_pdf_prior_plus_process_bu,l_prior_pdf,'Process BU','L')
      #else:
      #  process_likelihood_bu_weights=process_likelihood_weight_l(sampled_l_uniform_l)
      #  l_process_bu_likelihood=likelihood(sampled_l_uniform_l,process_likelihood_bu_weights,'Process BU','L')

    # ----------------------------------------------------------------------------
    # historical update
    # ----------------------------------------------------------------------------

    if no_hist:
      print ('No hist likelihood')
      s_hist_likelihood=np.full(n_bins-1,0.0)
      l_hist_likelihood=np.full(n_bins-1,0.0)
      total_hist_erf_prior_sample = np.full(nsamp,0.0)
    else:
      #sample size
      S=np.copy(sampled_s_prior)

      F2x = np.copy(f2co2_proc)

      l = -1*F2x/S

      #pattern effect parameter

      # updated value 0.6 Wm-2K-1, with a range from +0.3 Wm-2K-1 to +1.0 Wm-2K-1
      # (approximated as ±0.4 Wm-2K-1, 5-95% range)

      if no_pattern_effect:
        lam_p = np.zeros(nsamp)
      else:
        if cmip5_pattern_effect:
          lam_p = rnorm(nsamp,mean=-0.2,sd=0.4/1.64) # updated 27/3/20
        else:
          # Main version based on Andrews et al but with reduced lower bound
          # lam_p = rnorm(nsamp,mean=-0.6,sd=0.4/1.64)
          # updated to 5-95% range of 0.5+/-0.5 as requested by Kyle on 2/7/19
          lam_p_mu = -0.5
          lam_p_sd = 0.5/1.64
          if half_pattern_effect:
            # sensitivity test with halved uncertainty in lam_p
            lam_p = rnorm(nsamp,mean=lam_p_mu,sd=lam_p_sd/2.0)
          else:
            if fat_tails:
              lam_p = rfat(nsamp,mean=lam_p_mu,sd=lam_p_sd)
            else:
              # baseline calculation
              lam_p = rnorm(nsamp,mean=lam_p_mu,sd=lam_p_sd)

      # Use James Annan's pattern effect dependence
      if pattern_effect_dependence2:
        lam_p=np.copy(lam_p_h)

      #parameters for the historical calculation

      #histforc is a separate piece of code that returns the mean historical co2
      #forcing (co2_hist_mean) and an ensemble of other_forcing samples (F_other)

      (F_other,co2_hist_mean,total_aerosol_erf) = hist_forc(n_samples,
                              df_recent,df_baseline,uncs_neg,uncs_pos)

      #net imbalance = ocean heat uptake basically, figures from the text

      hist_imbalance = rnorm(nsamp,mean=0.6,sd=0.3/1.64)  # updated 19/2/19, checked 13/3/20

      co2_hist=co2_hist_mean*F2x/np.mean(F2x)

      total_hist_erf_prior_sample = F_other + co2_hist

      f_bin_boundaries=np.linspace(-10,10,n_bins+1)
      f_bin_width=f_bin_boundaries[1]-f_bin_boundaries[0]
      f_bin_centres=f_bin_boundaries[0:n_bins]+f_bin_width/2

      f=np.copy(co2_hist)
      flabel='co2_hist'
      f_pdf, f_bin_boundaries = np.histogram(f, bins=f_bin_boundaries)
      #f_pdf=normalise_pdf(f_pdf,f_bin_width)
      ([p_f_lt_1p5,p_f_gt_4p5,f5,f95,f_mean,f_mode])=ecs_cdf(f_pdf,f_bin_boundaries,f_bin_width,flabel)
      print(flabel,"=",f_mean,f5,f95)

      f=np.copy(total_aerosol_erf)
      flabel='total_aerosol_erf'
      f_pdf, f_bin_boundaries = np.histogram(f, bins=f_bin_boundaries)
      #f_pdf=normalise_pdf(f_pdf,f_bin_width)
      ([p_f_lt_1p5,p_f_gt_4p5,f5,f95,f_mean,f_mode])=ecs_cdf(f_pdf,f_bin_boundaries,f_bin_width,flabel)
      print(flabel,"=",f_mean,f5,f95)
      total_aerosol_erf_pdf = np.copy(f_pdf)

      f=np.copy(F_other)
      flabel='F_other'
      f_pdf, f_bin_boundaries = np.histogram(f, bins=f_bin_boundaries)
      #f_pdf=normalise_pdf(f_pdf,f_bin_width)
      ([p_f_lt_1p5,p_f_gt_4p5,f5,f95,f_mean,f_mode])=ecs_cdf(f_pdf,f_bin_boundaries,f_bin_width,flabel)
      print(flabel,"=",f_mean,f5,f95)

      f=np.copy(total_hist_erf_prior_sample)
      flabel='tot_forc'
      f_pdf, f_bin_boundaries = np.histogram(f, bins=f_bin_boundaries)
      ([p_f_lt_1p5,p_f_gt_4p5,f5,f95,f_mean,f_mode])=ecs_cdf(f_pdf,f_bin_boundaries,f_bin_width,flabel)
      print(flabel,"=",f_mean,f5,f95)

      if classic_hist_pdf:

        (s,l) = historical_classic_sample(F2x,co2_hist_mean,F_other,hist_imbalance,lam_p)

        s_pdf = sample2histogram(s, bin_boundaries, weights=equal_weights)
        l_pdf = sample2histogram(l, bin_boundaries, weights=equal_weights)

        s_pdf_updated = np.copy(s_pdf)
        l_pdf_updated = np.copy(l_pdf)

        s_hist_likelihood=np.full(n_bins,0.0)
        l_hist_likelihood=np.full(n_bins,0.0)

      else:
        #likelihood based on observed temperature rise of 0.91 with uncertainty of 0.15 at 90% = 1.64sd of gaussian, again taken from the text

        # **** Histrorical Likelihood ****

        # Calculate likelihood for S

        wt_hist = historical_likelihood(sampled_l_uniform_s,F2x,co2_hist_mean,F_other,hist_imbalance,lam_p)
        s_hist_likelihood=likelihood(sampled_s_uniform_s,wt_hist,'Historical','S')

        # Calculate likelihood for l

        wt_hist = historical_likelihood(sampled_l_uniform_l,F2x,co2_hist_mean,F_other,hist_imbalance,lam_p)
        l_hist_likelihood=likelihood(sampled_l_uniform_l,wt_hist,'Historical','L')

        # Calculate likelihood weights for Bayesian update

        wt_hist = historical_likelihood(l,F2x,co2_hist_mean,F_other,hist_imbalance,lam_p)

        #wt_hist = wt_hist / np.sum (wt_hist) #sum to unity to make proper weights

        #print(c("effective ensemble size based on historical",1/sum(wt_hist^2)))

        #another calculation ignoring the pattern effect

        # wt_hist_nopat = historical_likelihood(lambda,F2x,co2_hist_mean,F_other,hist_imbalance,lam_p*0)

        # wt_hist_nopat - wt_hist_nopat / sum (wt_hist_nopat) #sum to unity to make proper weights

        # print(c("effective ensemble size based on historical",1/sum(wt_hist_nopat^2)))

        weights = weights * wt_hist

        s_pdf_updated = sample2histogram(sampled_s_prior, bin_boundaries, weights=weights)
        #,inspect_stats='s_pdf_updated')

        #s_pdf_updated = normalise_pdf(s_pdf_updated,bin_width)

        #s_hist_pseudo_likelihood = pseudo_likelihood(s_pdf_updated,s_pdf,'Historical','S')

        # calculate the hist likelihood relative to the prior to reduce noise
        s_pdf_prior_plus_hist = sample2histogram(sampled_s_prior, bin_boundaries, weights=wt_hist)
        s_hist_pseudo_likelihood = pseudo_likelihood(s_pdf_prior_plus_hist,s_prior_pdf,'Historical','S')

        s_pdf=np.copy(s_pdf_updated)

        l_pdf_updated , bin_boundaries = np.histogram(sampled_l_prior, weights=weights, bins=bin_boundaries)
        #l_pdf_updated = normalise_pdf(l_pdf_updated,bin_width)

        l_hist_pseudo_likelihood = pseudo_likelihood(l_pdf_updated,l_pdf,'Historical','L',smooth=False)

        l_pdf=np.copy(l_pdf_updated)

    # ----------------------------------------------------------------------------
    # paleo update
    # ----------------------------------------------------------------------------

    S=np.copy(sampled_s_prior)

    nsamp=n_samples
    #shared parameters
    #S = runif(nsamp,min=0,max=10)

    F2x = np.copy(f2co2_proc)

    # Comments from James Annan's code supplied 16th March 2020

    # small changes to forcing as of 30 Jan 2020

    #comment from Gavin Schmidt:
    #Some research: If PI is 1850 as in the CMIP runs, CO2=284ppm, CH4=808ppb, N2O is 273ppb. LGM values are 190, 375, 200. Using updated forcing formulae from Etmanin et al (2016), the forcings are: CO2: -2.16, CH4: -0.37, N2O: -0.27, total -2.80 W/m2. I think that we should include a 1.45 factor for CH4 to account for O3 and strat WV components. So total would be: 2.97 W/m2.

    #my note: Etminan number is 0.57 of a doubling as F2x is provided separately
    #also note we are using 284ppm in the mPWP

    #note further as of 16 Mar the GHG forcings are all scaled up by 5%
    #this is via F2x for CO2 which also accounts for the scaling factors in warm periods
    #and directly in terms of the numbers for the minor gases at the LGM which
    #becomes: CO2: -2.27, CH4: -0.39*1.45 = -0.57, N2O: -0.28, total -3.12 W/m2.
    #with the CO2 as 0.57 of a doubling and the rest is 0.85

    #Updated as advised by James Annan 5/7/19

    #parameters for the paleo cold calculation

    # forcing excluding CO2
    # Updated values from James Annan 17/3/20
    # the numbers are ice, CH4, N2O, veg and dust in that order.
    if fat_tails:
      F_prime = rfat(nsamp,mean=-3.2-0.39*1.45-0.28-1.1-1.0,sd=2.)
    else:
      F_prime = rnorm(nsamp,mean=-3.2-0.39*1.45-0.28-1.1-1.0,sd=2.)

    #nonlinearity coefficient

    #this is big in the new formulation, potentially large change in lambda at 5C
    #as implemented it's a mean shift of 0.5 with an uncertainty of -0.5 to 1.5
    #at 2sd

    if fat_tails:
      alpha = rfat(nsamp,mean=-0.1,sd=0.1)
    else:
      alpha = rnorm(nsamp,mean=-0.1,sd=0.1)

    alpha_zero = alpha * 0 -0.0001 #nonzero to make the quadratic work!

    #pars for paleo hot

    #ppm CO2
    if fat_tails:
      ppm = rfat(nsamp,mean=375,sd=25)
    else:
      ppm = rnorm(nsamp,mean=375,sd=25)

    #CH4 scale-up

    #Gavin's factor of 1.4 as in his comments (accounting also for N2O)

    ch4_fac=rnorm(nsamp,mean=0.4,sd=0.1)

    #ESS scale-up
    if fat_tails:
      scalefac = rfat(nsamp,mean=0.5,sd=0.25)
    else:
      scalefac = rnorm(nsamp,mean=0.5,sd=0.25)

    # plot paleo hot forcing

    plt.xlim([-1,21])
    plt.ylim([0,1])
    plt.xlabel('S (K)')

    # code to plot breakdown of paleo hot calculations
    if 0:
      f_sample=(F2x)
      f_pdf=sample2histogram(f_sample,bin_boundaries,equal_weights)
      plt.plot(bin_centres,f_pdf,'k',linewidth=2)

      f_sample=(F2x * (1+scalefac))
      f_pdf=sample2histogram(f_sample,bin_boundaries,equal_weights)
      plt.plot(bin_centres,f_pdf,'r',linewidth=2)

      f_sample=((1+ch4_fac) * F2x * (1+scalefac) )
      f_pdf=sample2histogram(f_sample,bin_boundaries,equal_weights)
      plt.plot(bin_centres,f_pdf,'g',linewidth=2)

      f_pdf=sample2histogram(f_sample,bin_boundaries,equal_weights)
      plt.plot(bin_centres,f_pdf,'b',linewidth=2)

      dt_pdf=sample2histogram(dt_sample,bin_boundaries,equal_weights)
      plt.plot(bin_centres,dt_pdf,'m',linewidth=2)

    # code to plot alternative 'classical' hot paleo PDF

    if 0:
      dt_sample=np.random.normal(loc=3, scale=1.0, size=nsamp)
      f_sample=(np.log(ppm/284)/np.log(2) * (1+ch4_fac) * F2x * (1+scalefac) )
      s_sample=dt_sample*F2x/f_sample
      s_tmp_pdf=sample2histogram(s_sample,bin_boundaries,equal_weights)
      plt.plot(bin_centres,s_tmp_pdf,'c',linewidth=2)

      plt.savefig(out_dir+'/'+calc_id+'4.png')
      plt.close('all')

    s_pdf_pre_paleo=np.copy(s_pdf)

    # ----------------------------------------------------------------------------
    # paleo cold update
    # ----------------------------------------------------------------------------

    if no_paleo_cold:
      print ('No paleo cold likelihood')
      s_paleo_cold_likelihood=np.full(n_bins-1,0.0)
      l_paleo_cold_likelihood=np.full(n_bins-1,0.0)
    else:

      # **** Paleo Cold Likelihood ****

      # Calculate likelihood for S

      wt_paleo_cold = paleo_cold_likelihood(sampled_l_uniform_s,F2x,transfer,F_prime,alpha)
      #wt_paleo_cold = wt_paleo_cold / np.sum(wt_paleo_cold) #sum to unity to make proper weights
      s_paleo_cold_likelihood=likelihood(sampled_s_uniform_s,wt_paleo_cold,'Paleo Cold','S')

      # Calculate likelihood for l

      wt_paleo_cold = paleo_cold_likelihood(sampled_l_uniform_l,F2x,transfer,F_prime,alpha)
      #wt_paleo_cold = wt_paleo_cold / np.sum(wt_paleo_cold) #sum to unity to make proper weights
      l_paleo_cold_likelihood=likelihood(sampled_l_uniform_l,wt_paleo_cold,'Paleo Cold','L')

      # Calculate likelihood weights for Bayesian update

      wt_paleo_cold = paleo_cold_likelihood(sampled_l_prior,F2x,transfer,F_prime,alpha)
      #wt_paleo_cold = wt_paleo_cold / np.sum(wt_paleo_cold) #sum to unity to make proper weights

      weights = weights * wt_paleo_cold

      s_pdf_updated , bin_boundaries = np.histogram(sampled_s_prior, weights=weights, bins=bin_boundaries)
      #s_pdf_updated = normalise_pdf(s_pdf_updated,bin_width)

      #s_paleo_cold_pseudo_likelihood = pseudo_likelihood(s_pdf_updated,s_pdf,'Paleo Cold','S')

      # calculate the paleo col likelihood relative to the prior to reduce noise
      s_pdf_prior_plus_cold = sample2histogram(sampled_s_prior, bin_boundaries, weights=wt_paleo_cold)
      s_paleo_cold_pseudo_likelihood = pseudo_likelihood(s_pdf_prior_plus_cold,s_prior_pdf,'Paleo Cold','S')

      s_pdf=np.copy(s_pdf_updated)

      l_pdf_updated , bin_boundaries = np.histogram(sampled_l_prior, weights=weights, bins=bin_boundaries)
      #l_pdf_updated = normalise_pdf(l_pdf_updated,bin_width)

      #l_paleo_cold_pseudo_likelihood = pseudo_likelihood(l_pdf_updated,l_pdf,'Paleo Cold','L')

      # calculate the paleo cold likelihood relative to the prior to reduce noise
      l_pdf_prior_plus_cold = sample2histogram(sampled_l_prior, bin_boundaries, weights=wt_paleo_cold)
      l_paleo_cold_pseudo_likelihood = pseudo_likelihood(l_pdf_prior_plus_cold,l_prior_pdf,'Paleo Cold','L')

      l_pdf=np.copy(l_pdf_updated)

    # ----------------------------------------------------------------------------
    # paleo hot update
    # ----------------------------------------------------------------------------

    if no_paleo_hot:
      print ('No paleo hot likelihood')
      s_paleo_hot_likelihood=np.full(n_bins-1,0.0)
      l_paleo_hot_likelihood=np.full(n_bins-1,0.0)
    else:

      # **** Paleo Hot Likelihood ****

      # Calculate likelihood for S

      wt_paleo_hot = paleo_hot_likelihood(sampled_l_uniform_s,F2x,transfer,ppm,ch4_fac,scalefac)
      #wt_paleo_hot = wt_paleo_hot / np.sum(wt_paleo_hot) #sum to unity to make proper weights
      s_paleo_hot_likelihood=likelihood(sampled_s_uniform_s,wt_paleo_hot,'Paleo Hot','S')

      # Calculate likelihood for l

      wt_paleo_hot = paleo_hot_likelihood(sampled_l_uniform_l,F2x,transfer,ppm,ch4_fac,scalefac)
      #wt_paleo_hot = wt_paleo_hot / np.sum(wt_paleo_hot) #sum to unity to make proper weights
      l_paleo_hot_likelihood=likelihood(sampled_l_uniform_l,wt_paleo_hot,'Paleo Hot','L')

      # Calculate likelihood weights for Bayesian update

      wt_paleo_hot = paleo_hot_likelihood(sampled_l_prior,F2x,transfer,ppm,ch4_fac,scalefac)
      #wt_paleo_hot = wt_paleo_hot / np.sum(wt_paleo_hot) #sum to unity to make proper weights

      weights = weights * wt_paleo_hot

      s_pdf_updated = sample2histogram(sampled_s_prior, bin_boundaries, weights=weights)

      #s_paleo_hot_pseudo_likelihood = pseudo_likelihood(s_pdf_updated,s_pdf,'Paleo Hot','S')

      # calculate the paleo hot likelihood relative to the prior to reduce noise
      s_pdf_prior_plus_hot = sample2histogram(sampled_s_prior, bin_boundaries, weights=wt_paleo_hot)
      s_paleo_hot_pseudo_likelihood = pseudo_likelihood(s_pdf_prior_plus_hot,s_prior_pdf,'Paleo Hot','S')

      l_pdf_updated , bin_boundaries = np.histogram(sampled_l_prior, weights=weights, bins=bin_boundaries)
      #l_pdf_updated = normalise_pdf(l_pdf_updated,bin_width)

      #l_paleo_hot_pseudo_likelihood = pseudo_likelihood(l_pdf_updated,l_pdf,'Paleo Hot','L')

      # calculate the paleo hot likelihood relative to the prior to reduce noise
      l_pdf_prior_plus_hot = sample2histogram(sampled_l_prior, bin_boundaries, weights=wt_paleo_hot)
      l_paleo_hot_pseudo_likelihood = pseudo_likelihood(l_pdf_prior_plus_hot,l_prior_pdf,'Paleo Hot','L')

      l_pdf=np.copy(l_pdf_updated)

      s_paleo_pseudo_likelihood = pseudo_likelihood(s_pdf_updated,s_pdf_pre_paleo,'Paleo all','S')

      s_pdf=np.copy(s_pdf_updated)

      l_pdf_updated , bin_boundaries = np.histogram(sampled_l_prior, weights=weights, bins=bin_boundaries)
      #l_pdf_updated = normalise_pdf(l_pdf_updated,bin_width)

      l_paleo_pseudo_likelihood = pseudo_likelihood(l_pdf_updated,l_pdf,'Paleo all','L')

      l_pdf=np.copy(l_pdf_updated)

    if plot_paleo_likelihood:

        wt_paleo_cold = paleo_cold_likelihood(sampled_l_uniform_l,F2x,transfer,F_prime,alpha)
        wt_paleo_hot = paleo_hot_likelihood(sampled_l_uniform_l,F2x,transfer,ppm,ch4_fac,scalefac)
        l_paleo_likelihood=likelihood(sampled_l_uniform_l,wt_paleo_cold*wt_paleo_hot,'Paleo','L')

        wt_paleo_cold = paleo_cold_likelihood(sampled_l_uniform_s,F2x,transfer,F_prime,alpha)
        wt_paleo_hot = paleo_hot_likelihood(sampled_l_uniform_s,F2x,transfer,ppm,ch4_fac,scalefac)
        s_paleo_likelihood=likelihood(sampled_s_uniform_s,wt_paleo_cold*wt_paleo_hot,'Paleo','S')

    # Emergent contraints

    if no_process_ec:
      print ('No process EC likelihood')
      s_process_ec_likelihood=np.full(n_bins-1,0.0)
      l_process_ec_likelihood=np.full(n_bins-1,0.0)
    else:

      # Calculate S likelihood

      process_likelihood_ec_weights=process_likelihood_weight_emergent_constriants_s(sampled_l_uniform_s,F2x,lambda_ec)
      s_process_ec_likelihood=likelihood(sampled_s_uniform_s,process_likelihood_ec_weights,'Process EC','S')

      # Calculate l likelihood

      process_likelihood_ec_weights=process_likelihood_weight_emergent_constriants_s(sampled_l_uniform_l,F2x,lambda_ec)
      l_process_ec_likelihood=likelihood(sampled_l_uniform_l,process_likelihood_ec_weights,'Process EC','L')

      # Calculate weights for Bayesian update

      process_likelihood_ec_weights=process_likelihood_weight_emergent_constriants_s(sampled_l_prior,F2x,lambda_ec)

      weights = weights * process_likelihood_ec_weights

      s_pdf_updated , bin_boundaries = np.histogram(sampled_s_prior, weights=weights, bins=bin_boundaries)
      #s_pdf_updated = normalise_pdf(s_pdf_updated,bin_width)

      s_process_ec_pseudo_likelihood = pseudo_likelihood(s_pdf_updated,s_pdf,'Process EC','S')

      s_pdf=np.copy(s_pdf_updated)

      l_pdf_updated , bin_boundaries = np.histogram(sampled_l_prior, weights=weights, bins=bin_boundaries)
      #l_pdf_updated = normalise_pdf(l_pdf_updated,bin_width)

      l_pdf=np.copy(l_pdf_updated)

      l_process_ec_pseudo_likelihood = pseudo_likelihood(l_pdf_updated,l_pdf,'Process EC','L',smooth=False)

    if plot_process_likelihood:

        process_likelihood_ec_weights=process_likelihood_weight_emergent_constriants_s(sampled_l_uniform_l,F2x,lambda_ec)
        process_likelihood_bu_weights=process_likelihood_weight_l(sampled_l_uniform_l)
        l_process_likelihood=likelihood(sampled_l_uniform_l,process_likelihood_bu_weights*process_likelihood_ec_weights,'Process','L')

        process_likelihood_ec_weights=process_likelihood_weight_emergent_constriants_s(sampled_l_uniform_s,F2x,lambda_ec)
        process_likelihood_bu_weights=process_likelihood_weight_l(sampled_l_uniform_s)
        s_process_likelihood=likelihood(sampled_s_uniform_s,process_likelihood_bu_weights*process_likelihood_ec_weights,'Process','S')

    #-----------------------------------------------------------------------------
    # Calculate posterior for ECS and other quantities
    #-----------------------------------------------------------------------------

    ecs_pdf , bin_boundaries = np.histogram(ecs_sample, weights=weights, bins=bin_boundaries)

    print ('total_hist_erf_prior_sample.shape=',total_hist_erf_prior_sample.shape)
    print ('weights.shape=',weights.shape)

    total_hist_erf_prior,      f_bin_boundaries = np.histogram(total_hist_erf_prior_sample, bins=f_bin_boundaries)
    total_hist_erf_posterior , f_bin_boundaries = np.histogram(total_hist_erf_prior_sample, weights=weights, bins=f_bin_boundaries)

    #-----------------------------------------------------------------------------
    # update lambda posteriors
    #-----------------------------------------------------------------------------

    l_posterior=np.copy(l_pdf)

    posterior=np.copy(s_pdf)

    #-----------------------------------------------------------------------------
    # update partial sums
    #-----------------------------------------------------------------------------

    if toploop_index == 0:
      transfer_unweighted_prior_pdf_sum=np.copy(transfer_unweighted_prior_pdf)
    else:
      transfer_unweighted_prior_pdf_sum=transfer_unweighted_prior_pdf_sum+transfer_unweighted_prior_pdf

    if toploop_index == 0:
      transfer_weighted_prior_pdf_sum=np.copy(transfer_weighted_prior_pdf)
    else:
      transfer_weighted_prior_pdf_sum=transfer_weighted_prior_pdf_sum+transfer_weighted_prior_pdf

    if toploop_index == 0:
      total_hist_erf_posterior_sum=np.copy(total_hist_erf_posterior)
    else:
      total_hist_erf_posterior_sum=total_hist_erf_posterior_sum+total_hist_erf_posterior

    if toploop_index == 0:
      total_hist_erf_prior_sum=np.copy(total_hist_erf_prior)
    else:
      total_hist_erf_prior_sum=total_hist_erf_prior_sum+total_hist_erf_prior

    if toploop_index == 0:
      ecs_pdf_sum=np.copy(ecs_pdf)
    else:
      ecs_pdf_sum=ecs_pdf_sum+ecs_pdf

    if plot_process_likelihood:
      if toploop_index == 0:
        l_process_likelihood_sum=np.copy(l_process_likelihood)
      else:
        l_process_likelihood_sum=l_process_likelihood_sum+l_process_likelihood
      if toploop_index == 0:
        s_process_likelihood_sum=np.copy(s_process_likelihood)
      else:
        s_process_likelihood_sum=s_process_likelihood_sum+s_process_likelihood

    if plot_paleo_likelihood:
      if toploop_index == 0:
        l_paleo_likelihood_sum=np.copy(l_paleo_likelihood)
      else:
        l_paleo_likelihood_sum=l_paleo_likelihood_sum+l_paleo_likelihood
      if toploop_index == 0:
        s_paleo_likelihood_sum=np.copy(s_paleo_likelihood)
      else:
        s_paleo_likelihood_sum=s_paleo_likelihood_sum+s_paleo_likelihood

    if toploop_index == 0:
      s_pdf_sum=np.copy(s_pdf)
    else:
      s_pdf_sum=s_pdf_sum+s_pdf

    if toploop_index == 0:
      s_prior_pdf_sum=np.copy(s_prior_pdf)
    else:
      s_prior_pdf_sum=s_prior_pdf_sum+s_prior_pdf

    if toploop_index == 0:
      full_l_prior_pdf_sum=np.copy(full_l_prior_pdf)
    else:
      full_l_prior_pdf_sum=full_l_prior_pdf_sum+full_l_prior_pdf

    if toploop_index == 0:
      l_process_bu_likelihood_sum=np.copy(l_process_bu_likelihood)
    else:
      l_process_bu_likelihood_sum=l_process_bu_likelihood_sum+l_process_bu_likelihood

    if toploop_index == 0:
      l_process_bu_likelihood_sum=np.copy(l_process_bu_likelihood)
    else:
      l_process_bu_likelihood_sum=l_process_bu_likelihood_sum+l_process_bu_likelihood

    if toploop_index == 0:
      l_process_ec_likelihood_sum=np.copy(l_process_ec_likelihood)
    else:
      l_process_ec_likelihood_sum=l_process_ec_likelihood_sum+l_process_ec_likelihood

    if toploop_index == 0:
      l_hist_likelihood_sum=np.copy(l_hist_likelihood)
    else:
      l_hist_likelihood_sum=l_hist_likelihood_sum+l_hist_likelihood

    if toploop_index == 0:
      l_paleo_cold_likelihood_sum=np.copy(l_paleo_cold_likelihood)
    else:
      l_paleo_cold_likelihood_sum=l_paleo_cold_likelihood_sum+l_paleo_cold_likelihood

    if toploop_index == 0:
      l_paleo_hot_likelihood_sum=np.copy(l_paleo_hot_likelihood)
    else:
      l_paleo_hot_likelihood_sum=l_paleo_hot_likelihood_sum+l_paleo_hot_likelihood

    if toploop_index == 0:
      l_prior_sum=np.copy(l_prior)
    else:
      l_prior_sum=l_prior_sum+l_prior

    if toploop_index == 0:
      l_posterior_sum=np.copy(l_posterior)
    else:
      l_posterior_sum=l_posterior_sum+l_posterior

    if toploop_index == 0:
      s_process_bu_likelihood_sum=np.copy(s_process_bu_likelihood)
    else:
      s_process_bu_likelihood_sum=s_process_bu_likelihood_sum+s_process_bu_likelihood

    if toploop_index == 0:
      s_process_ec_likelihood_sum=np.copy(s_process_ec_likelihood)
    else:
      s_process_ec_likelihood_sum=s_process_ec_likelihood_sum+s_process_ec_likelihood

    if toploop_index == 0:
      s_hist_likelihood_sum=np.copy(s_hist_likelihood)
    else:
      s_hist_likelihood_sum=s_hist_likelihood_sum+s_hist_likelihood

    if toploop_index == 0:
      s_paleo_cold_likelihood_sum=np.copy(s_paleo_cold_likelihood)
    else:
      s_paleo_cold_likelihood_sum=s_paleo_cold_likelihood_sum+s_paleo_cold_likelihood

    if toploop_index == 0:
      s_paleo_hot_likelihood_sum=np.copy(s_paleo_hot_likelihood)
    else:
      s_paleo_hot_likelihood_sum=s_paleo_hot_likelihood_sum+s_paleo_hot_likelihood

    if toploop_index == 0:
      full_s_prior_pdf_sum=np.copy(full_s_prior_pdf)
    else:
      full_s_prior_pdf_sum=full_s_prior_pdf_sum+full_s_prior_pdf

    if toploop_index == 0:
      posterior_sum=np.copy(posterior)
    else:
      posterior_sum=posterior_sum+posterior

    # plot prior, likelihoods and posterior

    posterior=s_pdf_sum/(toploop_index+1)

    s_prior_pdf=s_prior_pdf_sum/(toploop_index+1)

    full_l_prior_pdf=full_l_prior_pdf_sum/(toploop_index+1)

    l_process_bu_likelihood=l_process_bu_likelihood_sum/(toploop_index+1)

    l_process_ec_likelihood=l_process_ec_likelihood_sum/(toploop_index+1)

    l_hist_likelihood=l_hist_likelihood_sum/(toploop_index+1)

    l_paleo_cold_likelihood=l_paleo_cold_likelihood_sum/(toploop_index+1)

    l_paleo_hot_likelihood=l_paleo_hot_likelihood_sum/(toploop_index+1)

    l_prior=l_prior_sum/(toploop_index+1)

    l_posterior=l_posterior_sum/(toploop_index+1)

    if plot_process_likelihood:
      s_process_likelihood=s_process_likelihood_sum/(toploop_index+1)
      l_process_likelihood=l_process_likelihood_sum/(toploop_index+1)

    if plot_paleo_likelihood:
      s_paleo_likelihood=s_paleo_likelihood_sum/(toploop_index+1)
      l_paleo_likelihood=l_paleo_likelihood_sum/(toploop_index+1)

    s_process_bu_likelihood=s_process_bu_likelihood_sum/(toploop_index+1)

    s_process_ec_likelihood=s_process_ec_likelihood_sum/(toploop_index+1)

    s_hist_likelihood=s_hist_likelihood_sum/(toploop_index+1)

    s_paleo_cold_likelihood=s_paleo_cold_likelihood_sum/(toploop_index+1)

    s_paleo_hot_likelihood=s_paleo_hot_likelihood_sum/(toploop_index+1)

    full_s_prior_pdf=full_s_prior_pdf_sum/(toploop_index+1)

    posterior=posterior_sum/(toploop_index+1)

    ecs_pdf=ecs_pdf_sum/(toploop_index+1)
    total_hist_erf_prior=total_hist_erf_prior_sum/(toploop_index+1)
    total_hist_erf_posterior=total_hist_erf_posterior_sum/(toploop_index+1)
    transfer_weighted_prior_pdf=transfer_weighted_prior_pdf_sum/(toploop_index+1)
    transfer_unweighted_prior_pdf=transfer_unweighted_prior_pdf_sum/(toploop_index+1)

    # normalise likelihoods

    if plot_process_likelihood:
      l_process_likelihood = normalise_likelihood(bin_centres,l_process_likelihood)
      s_process_likelihood = normalise_likelihood(bin_centres,s_process_likelihood)

    if plot_paleo_likelihood:
      l_paleo_likelihood = normalise_likelihood(bin_centres,l_paleo_likelihood)
      s_paleo_likelihood = normalise_likelihood(bin_centres,s_paleo_likelihood)

    l_process_bu_likelihood = normalise_likelihood(bin_centres,l_process_bu_likelihood)
    l_process_ec_likelihood = normalise_likelihood(bin_centres,l_process_ec_likelihood)
    l_hist_likelihood = normalise_likelihood(bin_centres,l_hist_likelihood)
    l_paleo_cold_likelihood = normalise_likelihood(bin_centres,l_paleo_cold_likelihood)
    l_paleo_hot_likelihood = normalise_likelihood(bin_centres,l_paleo_hot_likelihood)

    s_process_bu_likelihood = normalise_likelihood(bin_centres,s_process_bu_likelihood)
    s_process_ec_likelihood = normalise_likelihood(bin_centres,s_process_ec_likelihood)
    s_hist_likelihood = normalise_likelihood(bin_centres,s_hist_likelihood)
    s_paleo_cold_likelihood = normalise_likelihood(bin_centres,s_paleo_cold_likelihood)
    s_paleo_hot_likelihood = normalise_likelihood(bin_centres,s_paleo_hot_likelihood)

    if toploop_index % plot_frequency == 0:

      dumpfile=out_dir+'/'+calc_id+'.means.'+str(toploop_index).zfill(6)+'.joblib'
      print ('dumpfile=',dumpfile)
      # Saving the objects:
      dump([transfer_unweighted_prior_pdf,transfer_weighted_prior_pdf,
            total_hist_erf_posterior,total_hist_erf_prior,ecs_pdf, posterior,
            n_bins, bin_boundaries, bin_centres, bin_width, n_samples, s_pdf,
            s_prior_pdf, full_l_prior_pdf, l_process_bu_likelihood,
            l_process_ec_likelihood, l_hist_likelihood, l_paleo_cold_likelihood,
            l_paleo_hot_likelihood, l_prior, l_posterior, s_process_bu_likelihood,
            s_process_ec_likelihood, s_hist_likelihood, s_paleo_cold_likelihood,
            s_paleo_hot_likelihood, full_s_prior_pdf], dumpfile)

      dumpfile=out_dir+'/'+calc_id+'.sums.'+str(toploop_index).zfill(6)+'.joblib'
      print ('dumpfile=',dumpfile)
      # Saving the objects:
      dump([toploop_index,transfer_unweighted_prior_pdf_sum,transfer_weighted_prior_pdf_sum,
            total_hist_erf_posterior_sum,total_hist_erf_prior_sum,ecs_pdf_sum,s_pdf_sum,
            s_prior_pdf_sum, full_l_prior_pdf_sum, l_process_bu_likelihood_sum,
            l_process_ec_likelihood_sum, l_hist_likelihood_sum, l_paleo_cold_likelihood_sum,
            l_paleo_hot_likelihood_sum, l_prior_sum, l_posterior_sum, s_process_bu_likelihood_sum,
            s_process_ec_likelihood_sum, s_hist_likelihood_sum, s_paleo_cold_likelihood_sum,
            s_paleo_hot_likelihood_sum, full_s_prior_pdf_sum] , dumpfile)

      # Getting back the objects:
      #with open('objs.pkl') as f:  # Python 3: open(..., 'rb')
      # obj0, obj1, obj2 = pickle.load(f)

      print('XX posterior',toploop_index)

      print('XX start of the end',toploop_index)

      biplot=0

      # plot prior, likelihoods and posterior

      import matplotlib.pyplot as plt

      posterior=s_pdf_sum/(toploop_index+1)

      plt.rcParams['figure.figsize'] = [10, 7]

      plt.close('all')

      if biplot:
        # Two subplots, the axes array is 1-d
        f, axarr = plt.subplots(1,2, sharey=True)

      color='k'

      # **** lambda plot ******

      l_prior=np.copy(full_l_prior_pdf)
      l_prior_label='Lambda '+prior_label

      xpos=-11
      pos=0.8
      dy=0.07

      color='k'
      #plt.plot([1.5,1.5],yrange,color=color)
      #plt.plot([4.5,4.5],yrange,color=color)

      #pos=2
      #color='k'
      #label='s_prior_pdf'
      #plt.text(-1,1-0.05*(pos+1),label,color='r')
      color='r'
      #plt.plot(bin_centres,s_prior_pdf,color)

      yrange=[0,2]
      if not plot_posterior:
        yrange=[0,1]

      if biplot:
        plt=axarr[0]
        plt.set_xlim([-6,1])
        plt.set_ylim(yrange)
        plt.set_xlabel('L (Wm-2K-1)')
      else:
        plt.xlim([-6,1])
        plt.ylim(yrange)
        plt.xlabel('L (Wm-2K-1)')

      if plot_prior:
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),l_prior_label,color=color)
        #plt.plot(bin_centres,l_prior_pdf,color,linestyle='dotted')
        plt.plot(bin_centres,normalise_pdf(full_l_prior_pdf,bin_width),color)

      pos=3

      if not no_process_bu:
        pos=pos+1
        color='g'
        label='BU Process Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,l_process_bu_likelihood,color,linewidth=2)

        # also plot total process likleihood explicitly

        #(lambda_mu,lambda_sd)=process_lambda_parameters()

        #process_lambda_sample = np.random.normal(loc=lambda_mu, scale=lambda_sd, size=nsamp*10)
        #process_lambda_likelihood = sample2histogram(process_lambda_sample,bin_boundaries,np.full(nsamp*10,1.0))
        #norm_process_lambda_likelihood = normalise_likelihood(bin_centres,process_lambda_likelihood)
        #
        #plt.plot(bin_centres,norm_process_lambda_likelihood,'g',linestyle='dashed')
        #plt.plot(bin_centres,dnorm(bin_centres, lambda_mu,lambda_sd),'g',linestyle='dotted')

      if not no_process_ec:
        pos=pos+1
        color='y'
        label='EC Process Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,l_process_ec_likelihood,color,linewidth=2)

      if not no_hist:
        pos=pos+1
        color='b'
        label='Historical Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,l_hist_likelihood,color,linewidth=2)

      if not no_paleo_cold:
        pos=pos+1
        color='c'
        label='Cold Paleo Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,l_paleo_cold_likelihood,color,linewidth=2)

      if not no_paleo_hot:
        pos=pos+1
        color='m'
        label='Hot Paleo Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,l_paleo_hot_likelihood,color,linewidth=2)

      if plot_process_likelihood:
        pos=pos+1
        color='b'
        label='Process Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,l_process_likelihood,color,linewidth=2)

      if plot_paleo_likelihood:
        pos=pos+1
        color='b'
        label='Paleo Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,l_paleo_likelihood,color,linewidth=2)

      color='r'
      print ('statistics for l_prior:')

      if plot_prior:
        #plt.plot(bin_centres,l_prior,linewidth=2,color=color)

        ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(full_l_prior_pdf,bin_boundaries,bin_width,'full_l_prior_pdf')

        #plt.text(xpos,yrange[1]*0.85,'P(ECS < 1.5) = '+"{:.2f}".format(p_ecs_lt_1p5/100.0),color=color)
        #plt.text(xpos,yrange[1]*0.8,'P(ECS > 4.5) = '+"{:.2f}".format(p_ecs_gt_4p5/100.0),color=color)
        plt.text(xpos,yrange[1]*0.85,'L 5th %ile = '+"{:.1f}".format(ecs5),color=color)
        plt.text(xpos,yrange[1]*0.8,'L 95th %ile = '+"{:.1f}".format(ecs95),color=color)

      print ('statistics for l_posterior:')

      if plot_posterior:
        color='k'
        plt.plot(bin_centres,normalise_pdf(l_posterior,bin_width),linewidth=2,color=color)
        plt.text(xpos,yrange[1]*0.5,l_posterior_label,color=color)

        ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(l_posterior,bin_boundaries,bin_width,'L Posterior')

        #plt.text(xpos,yrange[1]*0.45,'P(ECS < 1.5) = '+"{:.2f}".format(p_ecs_lt_1p5/100.0),color=color)
        #plt.text(xpos,yrange[1]*0.4,'P(ECS > 4.5) = '+"{:.2f}".format(p_ecs_gt_4p5/100.0),color=color)
        plt.text(xpos,yrange[1]*0.45,'L 5th %ile = '+"{:.1f}".format(ecs5),color=color)
        plt.text(xpos,yrange[1]*0.4,'L 95th %ile = '+"{:.1f}".format(ecs95),color=color)

      print ('xx l posterior',toploop_index)
      #if biplot == 0:
      plt.savefig(out_dir+'/'+calc_id+'l_stats.png')
      plt.close('all')

      # **** S plot ******

      yrange=[0,1]
      if biplot:
        plt=axarr[1]
        plt.set_xlim([-1,21])
        plt.set_ylim(yrange)
        plt.set_xlabel('S (K)')
      else:
        plt.xlim([-1,21])
        plt.ylim(yrange)
        plt.xlabel('S (K)')

      xpos=8
      pos=3

      if not no_process_bu:
        pos=pos+1
        color='g'
        label='BU Process Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,s_process_bu_likelihood,color,linewidth=2)

      if not no_process_ec:
        pos=pos+1
        color='y'
        label='EC Process Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,s_process_ec_likelihood,color,linewidth=2)


      if not no_hist:
        pos=pos+1
        color='b'
        label='Historical Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,s_hist_likelihood,color,linewidth=2)

      if not no_paleo_cold:
        pos=pos+1
        color='c'
        label='Cold Paleo Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,s_paleo_cold_likelihood,color,linewidth=2)

      if not no_paleo_hot:
        pos=pos+1
        color='m'
        label='Hot Paleo Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,s_paleo_hot_likelihood,color,linewidth=2)

      if plot_process_likelihood:
        pos=pos+1
        color='b'
        label='Process Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,s_process_likelihood,color,linewidth=2)

      if plot_paleo_likelihood:
        pos=pos+1
        color='b'
        label='Paleo Likelihood'
        plt.text(xpos,yrange[1]*(1-0.05*(pos+1)),label,color=color)
        plt.plot(bin_centres,s_paleo_likelihood,color,linewidth=2)

      color='r'
      print ('statistics for prior:')

      if plot_prior:
        #plt.plot(bin_centres,s_prior_pdf,linewidth=2,color=color,linestyle='dotted')
        plt.plot(bin_centres,normalise_pdf(full_s_prior_pdf,bin_width),color)

        plt.text(xpos,yrange[1]*0.9,prior_label,color=color)

        ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(s_prior_pdf,bin_boundaries,bin_width,'s_prior_pdf')

        ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(full_s_prior_pdf,bin_boundaries,bin_width,'full_s_prior_pdf')

        #plt.text(xpos,yrange[1]*0.85,'P(ECS < 1.5) = '+"{:.2f}".format(p_ecs_lt_1p5/100.0),color=color)
        #plt.text(xpos,yrange[1]*0.8,'P(ECS > 4.5) = '+"{:.2f}".format(p_ecs_gt_4p5/100.0),color=color)
        plt.text(xpos,yrange[1]*0.85,'ECS 5th %ile = '+"{:.1f}".format(ecs5),color=color)
        plt.text(xpos,yrange[1]*0.8,'ECS 95th %ile = '+"{:.1f}".format(ecs95),color=color)

      print ('statistics for posterior:')

      if plot_posterior:
        color='k'
        # for this working plot we call normalise_likelihood which scales pdf to peak value as can be spiky at the start
        plt.plot(bin_centres,normalise_likelihood(bin_centres,posterior) ,linewidth=2,color=color)
        plt.text(xpos,yrange[1]*0.5,s_posterior_label,color=color)

        ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(posterior,bin_boundaries,bin_width,'Loopsum Posterior')
        ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(ecs_pdf,bin_boundaries,bin_width,'ecs_pdf')
        ([p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs95,ecs_mean,ecs_mode])=ecs_cdf(total_hist_erf_posterior_sum,f_bin_boundaries,f_bin_width,'total_hist_erf_posterior_sum')

        #plt.text(xpos,yrange[1]*0.45,'P(ECS < 1.5) = '+"{:.2f}".format(p_ecs_lt_1p5/100.0),color=color)
        #plt.text(xpos,yrange[1]*0.4,'P(ECS > 4.5) = '+"{:.2f}".format(p_ecs_gt_4p5/100.0),color=color)
        plt.text(xpos,yrange[1]*0.45,'ECS 5th %ile = '+"{:.1f}".format(ecs5),color=color)
        plt.text(xpos,yrange[1]*0.4,'ECS 95th %ile = '+"{:.1f}".format(ecs95),color=color)

        print ('xx s posterior')

        if calc_id == 'US2':
          plt.plot(bin_centres,normalise_pdf(full_importance_density_scaling,bin_width) ,linewidth=2,color='g')

      #if biplot == 0:
      plt.savefig(out_dir+'/'+calc_id+'s_stats.png')
      plt.close('all')

      if calc_id == 'PLOT_ERF' or calc_id == 'PLOT_ERF_AR5' or calc_id == 'PLOT_ERF_BELLCON':
        plt.rcParams['figure.figsize'] = [15, 6]
        plt.rcParams.update({'font.size': 15})
        plt.close('all')
        yrange=[0,1]
        xrange=[-1,4]
        plt.xlim(xrange)
        plt.ylim(yrange)
        plt.xlabel('Total ERF $(Wm^{-2}$)')
        plt.ylabel('PDF for Total ERF $(W^{-1}m^{2}$)')
        #plt.xticks(np.arange(1, 13, step=1))
        total_hist_erf_prior=normalise_pdf(total_hist_erf_prior,f_bin_width)
        plt.plot(f_bin_centres,total_hist_erf_prior,'b',linewidth=4,label='Total ERF')

        legend = plt.legend(loc='upper right', fontsize='medium')
        plt.savefig(out_dir+'/'+'fig4.1_debug.png')

      print ('XX posterior')

      print('XX end')
