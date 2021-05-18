import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import os

def bellouin_ringberg_stats_fix_200311_ERFaer(nsamples,
                        doplots=False, fig_dir='./figs/'):
  # Original code provided by # Ed Gryspeerdt (Imperial College London)

  # Modified by Nicolas Bellouin (University of Reading)
  # to reproduce forcings from # published version of Bellouin et al,
  # 2020:
  # Bounding global aerosol radiative forcing of climate change
  # [1]https://doi.org/10.1029/2019RG000660

  # Minor modifications by Mark Webb, Met Office,
  # for use with WCRP ECS Assessment code

  # Minor modifications by Qun Liu, Univ. of Exeter,
  # for reproducing the Fig. 5 of Myers et al. (2021)

  # Distribution type
  # u - uniform, but arguments are 16th and 84th percentile
  #       assumes uniform outside these values (no very large values)
  # v - uniform, next two arguments are upper and lower bounds
  # c - constant value
  # n - normal distribution with bounds as 2 sigma width around an average mean

  d = {}
  d['rsdt'] =          ['c', -340]

  # Aerosol variables
  d['dtau'] =          ['u', 0.02, 0.04]
  d['tau'] =           ['u', 0.13, 0.17]

  # RFari terms
  d['S_tau'] =         ['u', -20,  -27]
  d['RFari_cloudy'] =  ['u', -0.1, 0.1]

  # RAari
  d['dR_dRatm'] =      ['u', -0.1, -0.3]
  d['dRatm_dtau'] =    ['u', 17,   35  ]

  # Cloud fractions
  d['c_tau'] =         ['u', 0.59, 0.71]
  d['c_N'] =           ['u', 0.19, 0.29]
  d['c_L'] =           ['u', 0.21, 0.29]
  d['c_C'] =           ['u', 0.59, 1.07]

  # Nd sensitivity to aerosol, and adjustment terms
  d['beta_N_tau'] =    ['u',  0.3,   0.8]
  d['beta_L_N'] =      ['u', -0.36, -0.011]
  d['beta_C_N'] =      ['u',  0,     0.1]

  # Cloud albedo terms
  d['S_N'] =           ['u', -26,  -27]
  d['S_L'] =           ['u', -54,  -56]
  d['S_C'] =           ['u', -91, -153]

  ############################################
  # Create the distributions and store in ds #
  ############################################

  ds = {}
  for name in d.keys():
      if d[name][0] == 'c':
          ds[name] = d[name][1]
      elif d[name][0] == 'u':
          diff = 16 * (d[name][2] - d[name][1]) / (84 - 16)
          ds[name] = np.random.uniform(d[name][1]-diff, d[name][2]+diff, nsamples)
      elif d[name][0] == 'v':
          ds[name] = np.random.uniform(d[name][1], d[name][2], nsamples)
      elif d[name][0] == 'n':
          cent = (d[name][2] + d[name][1]) / 2
          spread = (d[name][2] - d[name][1]) / 2
          ds[name] = np.random.uniform(cent, spread, nsamples)

  ############################
  # The Ringberg Equation!!! #
  ############################
  rfari = ds['dtau'] * ds['S_tau'] * (1 - ds['c_tau']) + ds['RFari_cloudy']
  rfari_adj = ds['dtau'] * ds['dR_dRatm'] * ds['dRatm_dtau']

  ds['dlntau'] = ds['dtau'] / ds['tau']
  deltan = ds['dlntau'] * ds['beta_N_tau']

  rfaci = ds['dlntau'] * ds['beta_N_tau'] * ds['S_N'] * ds['c_N']
  erfaci_L = ds['dlntau'] * ds['beta_N_tau'] * ds['beta_L_N'] * ds['S_L'] * ds['c_L']
  erfaci_C = ds['dlntau'] * ds['beta_N_tau'] * ds['beta_C_N'] * ds['S_C'] * ds['c_C']

  ERFaer = rfari + rfari_adj + rfaci + erfaci_L + erfaci_C

  #########################
  # Output and formatting #
  #########################

  output = [('dtau/tau', ds['dlntau']),
            ('deltaN', deltan),
            ('RFari', rfari),
            ('RFari_adj', rfari_adj),
            ('ERFari', rfari + rfari_adj),
            (None),
            ('RFaci', rfaci),
            ('ERFaci_L', erfaci_L),
            ('ERFaci_C', erfaci_C),
            ('ERFaci', rfaci + erfaci_L + erfaci_C),
            (None),
            ('ERFaer', ERFaer)]
            #('ERFaer', rfari + rfari_adj + rfaci + erfaci_L + erfaci_C)]

  print('')
  print('{: >10}  {: >6} {: >6} {: >6}'.format('', 'Mean', 'Lower', 'Upper'))
  for op in output:
      if op:
          print('{: >10}: {:6.2f} {:6.2f} {:6.2f}'.format(
              op[0],
              np.mean(op[1]),
              np.percentile(op[1], [16])[0],
              np.percentile(op[1], [84])[0]))
      else:
          print('')

  #doplots = #False
  if doplots:
    # note that the hard coded range bars in the plots below are not up to
    # date with the final numbers in Bellouin et al 2020
    # These aren't used in the WCRP ECS assessment

    bin_width = 0.01

    fig, axs = plt.subplots(2, 2) #, figsize=(11.69,8.27))

    # ARI
    axs[0,0].hist(rfari, density=True, #histtype='step', # normed=True,
             bins=np.arange(-4.5, 1.5, bin_width), linewidth=2, color='b',
             linestyle='--')
    axs[0,0].hist(rfari + rfari_adj, density=True, #normed=True,histtype='step',
             bins=np.arange(-4.5, 1.5, bin_width),linewidth=2, color='b')
    axs[0,0].set_title('Aerosol-radiation interactions')
    axs[0,0].set_xlabel(r'Radiative Forcing (Wm$^{-2}$)')
    axs[0,0].set_xlim(-4, 0.5)
    axs[0,0].errorbar(-0.45, 5.0, xerr=0.25, linewidth=2,capsize=5,color='b')
    axs[0,0].errorbar(-0.45, 5.5, xerr=0.5, fmt='o',linewidth=2,capsize=5,color='k')
    axs[0,0].grid(True)

    # ACI
    axs[0,1].hist(rfaci, density=True,  #histtype='step', #normed=True,
             bins=np.arange(-4.5, 1.5, bin_width),linewidth=2,color='b',
             linestyle='--')
    axs[0,1].hist(rfaci+erfaci_L+erfaci_C, density=True,  #histtype='step', #normed=True,
             bins=np.arange(-4.5, 1.5, bin_width),linewidth=2,color='b')
    axs[0,1].set_title('Aerosol-cloud interactions')
    axs[0,1].set_xlabel(r'Radiative Forcing (Wm$^{-2}$)')
    axs[0,1].set_xlim(-4, 0.5)
    axs[0,1].errorbar(-3.10, 1.2, xerr=[[0.0],[3.0]], linewidth=2,capsize=5,color='b')
    axs[0,1].errorbar(-0.45, 1.3, xerr=[[0.75],[0.45]], fmt='o',linewidth=2,capsize=5,color='k')
    axs[0,1].grid(True)

    # Total
    axs[1,0].hist(rfari + rfaci,  density=True, #histtype='step', #normed=True,
             bins=np.arange(-4.5, 1.5, bin_width),linewidth=2,color='b',
             linestyle='--')
    axs[1,0].hist(rfari + rfari_adj + rfaci+erfaci_L+erfaci_C,  density=True, #histtype='step', #normed=True,
             bins=np.arange(-4.5, 1.5, bin_width),linewidth=2,color='b')
    axs[1,0].set_title('Total aerosol')
    axs[1,0].set_xlabel(r'Radiative Forcing (Wm$^{-2}$)')
    axs[1,0].set_xlim(-4, 0.5)
    axs[1,0].errorbar(-3.6, 1.2, xerr=[[0.0],[3.2]], linewidth=2,capsize=5,color='b')
    axs[1,0].errorbar(-0.9, 1.3, xerr=[[1.0],[0.8]], fmt='o',linewidth=2,capsize=5,color='k')

    # IPCC PDF, see Supplementary material of Chapter 8, page 8SM-11
    alpha = 0.2425
    f = 1.645
    p05 = -1.9
    p95 = -0.1
    best = -0.9
    x0 = (((4./5.)*math.sqrt(2.0/math.pi) * (p05+p95))-((2*alpha+f)*best)) / (((8./5.)*math.sqrt(2.0/math.pi))-(2*alpha+f))
    sigmap = ((alpha+f) * (p95-x0)+alpha * (x0-p05)) / ((alpha+f) * (alpha+f)-alpha*alpha)
    sigmam = ((alpha+f) * (x0-p05)+alpha * (p95-x0)) / ((alpha+f) * (alpha+f)-alpha*alpha)
    xm = np.arange(-4.5, x0, bin_width)
    ipccm=1/math.sqrt(2*math.pi)*2.0/(sigmap+sigmam)*np.exp(-((xm-x0) * (xm-x0)) / (2*sigmam*sigmam))
    axs[1,0].plot(xm,ipccm,linewidth=3,color='k')
    xp=np.arange(x0,1.5,bin_width)
    ipccp=1/math.sqrt(2*math.pi)*2.0/(sigmap+sigmam)*np.exp(-((xp-x0) * (xp-x0)) / (2*sigmap*sigmap))
    axs[1,0].plot(xp,ipccp,linewidth=3,color='k')

    axs[1,0].grid(True)
    rect = patches.Rectangle((-4,0),2,1.4, color='r',alpha=0.75,zorder=0)
    axs[1,0].add_patch(rect)
    rect = patches.Rectangle((-4,0),2.4,1.4, color='r',alpha=0.5,zorder=0)
    axs[1,0].add_patch(rect)
    rect = patches.Rectangle((0,0),0.5,1.4, color='gold',alpha=0.75,zorder=0)
    axs[1,0].add_patch(rect)

    # Empty panel
    axs[1,1].axis('off')

    #fig.subplots_adjust(hspace=0.3)
    fig.tight_layout()

    plt.savefig(os.path.join(fig_dir, "ringberg_stats.pdf"), dpi=200)
    plt.show()

  return (ERFaer)


if __name__ == '__main__':
    fig_dir = './figs'
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    nsamples = 10000
    plt.close()
    ERFaer = bellouin_ringberg_stats_fix_200311_ERFaer(nsamples,
             doplots=True, fig_dir=fig_dir)
