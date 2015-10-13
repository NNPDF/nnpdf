# nnpydf - nh 08/14

import commondata
import os, glob, re

def markdown_data(data, prefix):
  print "Exporting Markdown for: " + data.setname

  # Fetch available plots
  plots = [f for f in os.listdir(prefix)]

  targetfile = prefix + "index.md"
  file = open(targetfile,'w')

  # Header
  file.write(data.setname+'\n')
  file.write("="*len(data.setname)+'\n')

  # noKin_plot
  file.write('#### Plot vs Datapoint \n')
  file.write('[!['+data.setname+' datapoints'+ ']('+data.setname+'.png)]('+data.setname+'.pdf) \n\n')
  
  file.write('[Return to Index](../index.html)\n\n')
  file.write('------------- \n')
  
  # write total kin plots
  file.write('#### Plot vs Kinematics (collated bins) \n')
  file.write('###### n.b bins are scaled by a factor of 2^i where i is the bin index  \n')
  for f in plots:
    for i in range(0, 200): # Must be a better way than this!
      if f == (data.setname+'_'+str(i)+'.pdf'):
        plot = data.setname+'_'+str(i)+'.png'
        file.write('[!['+data.setname+'_'+str(i)+ ']('+plot+')]('+f+')\n')

  file.write('      \n')
  file.write('[Return to Index](../index.html)\n\n')
  file.write('------------- \n')

  # write total kin plots - ratios
  file.write('#### Ratio plot vs Kinematics (collated bins) \n')
  for f in plots:
    for i in range(0, 200): # Must be a better way than this!
      if f == (data.setname+'_'+str(i)+'_R.pdf'):
        plot = data.setname+'_'+str(i)+'_R.png'
        file.write('[!['+data.setname+'_'+str(i)+ ']('+plot+')]('+f+')\n')

  file.write('      \n')
  file.write('[Return to Index](../index.html)\n\n')
  file.write('------------- \n')

  file.write('#### Plot vs Kinematics (individual bins) \n')
  for f in plots:
    for i in range(0, 200):
      for j in range(0, 200):
        if f == (data.setname+'_'+str(j)+'_'+str(i)+'.pdf'):
          plot = data.setname+'_'+str(j)+'_'+str(i)+'.png'
          file.write('[!['+data.setname+'_'+str(j)+'_'+str(i)+ ']('+plot+')]('+f+')\n')

  file.write('      \n')
  file.write('[Return to Index](../index.html)\n\n')
  file.write('------------- \n')

  file.close()