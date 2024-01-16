import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

frame = int(input('Frame number : '))
fname = 'q1_'+str(frame).zfill(3)+'.r4'

#nx = int(input('  - nx = '))
#ny = int(input('  - ny = '))

n = int(input('  - n = '))
nx = n
ny = n

in_file = open(fname,'r')

array = np.fromfile(fname,dtype=np.float32)

a1 = array[1:nx*ny+1].reshape(nx,ny)
a2 = array[nx*ny+2:].reshape(nx,ny)

a1 = np.transpose(np.array(a1))
a2 = np.transpose(np.array(a2))

a1max = np.max(a1)
a1min = np.min(a1)
lev1 = max(abs(a1max),abs(a1min))
a2max = np.max(a2)
a2min = np.min(a2)
lev2 = max(abs(a2max),abs(a2min))

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ima = ax.imshow(a1,cmap='seismic',vmin=-lev1,vmax=lev1)
plt.colorbar(ima)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

fig2 = plt.figure()
bx = fig2.add_subplot(111)
imb = bx.imshow(a2,cmap='seismic',vmin=-lev2,vmax=lev2)
plt.colorbar(imb)
bx.set_xticklabels([])
bx.set_yticklabels([])
bx.xaxis.set_ticks_position('none')
bx.yaxis.set_ticks_position('none')

fn1 = 'im1_'+str(frame).zfill(3)+'.png'
fn2 = 'im2_'+str(frame).zfill(3)+'.png'

fig1.savefig(fn1,format='png',bbox_inches='tight')
fig2.savefig(fn2,format='png',bbox_inches='tight')

plt.show()    

