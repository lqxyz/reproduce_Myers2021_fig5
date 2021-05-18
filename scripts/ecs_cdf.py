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

import numpy as np
from normalise_pdf import *

def ecs_cdf(count_in,bin_boundaries,bin_width,label):

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

  print ( 'Prior;P(ECS <1);P(ECS <1.5);P(ECS <2);P(ECS >4);P(ECS >4.5);P(ECS >6);\
          P(1.5 -4.5);5th pile;10th pile;17th pile;20th pile;25th pile;50th pile;\
          75th pile;80th pile;83rd pile;90th pile;95th pile;Mode;Mean' )
  qlabel=label
  f='10.2f'
  f=''
  print ("qlabel=",qlabel, p_ecs_lt_1, p_ecs_lt_1p5, p_ecs_lt_2,100-p_ecs_lt_4,
          100-p_ecs_lt_4p5, 100-p_ecs_lt_6,p_ecs_lt_4p5-p_ecs_lt_1p5,ecs5,
          ecs10, ecs17, ecs20,ecs25,ecs50, ecs75,ecs80,ecs83, ecs90, ecs95,
          ecs_mode,ecs_mean)

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

  f='10.2f'
  print (format('median and 5-95% range for '+label+':','10s'),
         format(ecs50 ,f),format(ecs5 ,f),'-',format(ecs95 ,f))

  slabel2=slabel+'.txt'
  slabel2=slabel2.replace(' ','_')
  slabel2=slabel2.replace('-','_')
  slabel2=slabel2.replace('/','_')

  #file = open(slabel2,"w")
  #file.write(bin_boundaries)
  #file.write(count)
  #file.close()

  return (p_ecs_lt_1p5,p_ecs_gt_4p5,ecs5,ecs17,ecs83,ecs95,ecs_mean,ecs_mode)
