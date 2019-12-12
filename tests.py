import numpy as np

import ofu
from ofu_tools import angular_distance, ang_distance, doublet_average

# a = np.arange(10)
# t = np.array([0.,1000.,1020.,1070.,2000.,2001.,2005.,2090.,3000.])
#
# print( np.diff(t,n=1) )
# print( np.diff(t,n=2) )
#
# print( np.diff(t,n=3) )
#
# print(t[np.append(np.diff(t,n=1)<1000.,np.array([False]) )])
#
data = []
for year in range(2011,2019):
   data.append(np.load('../icecube_datasets/gfu/version-002-p04/IC86_{}_data.npy'.format(year)) )
data = np.concatenate(data)
#
# print(data.dtype)
# m=4
# for i in range(m,0,-1):
#     for j in range(i):
#         print(m-i,-i,' ',m-j,-j)
#
# a= [0,1,2]
# print(a[1:])

test = ofu.multiplets(data,100,3.5,floor=0.2)
print(test.floor)
print(test.n_multiplets)
print(test.get_doublets)
doublets = test.get_doublets
print(doublets['t1'])
print(np.sum(np.diff(doublets['t1'])<0) )
#
# multiplicity = 4
# space_diff = angular_distance(data['ra'][:-(multiplicity-1)],data['dec'][:-(multiplicity-1)],data['ra'][(multiplicity-1):],data['dec'][(multiplicity-1):])
# space_diff_2 = ang_distance(data[:-(multiplicity-1)],data[(multiplicity-1):])
#
# print(space_diff)
# print(space_diff_2)
# print(angular_distance(data['ra'][0],data['dec'][0],data['ra'][multiplicity-1],data['dec'][multiplicity-1]))
# print(angular_distance(data['ra'][1],data['dec'][1],data['ra'][multiplicity],data['dec'][multiplicity]))

# a = np.array([1,2,3,4])
# b = np.array([5,6,7,8])
# c=np.column_stack((a,b))
#
# print(c)
# print(np.mean(c,axis=0) )

# av = doublet_average(data[0:4],data[1:5],theta_a=0.9)
# print(av)
# print(av.dtype)
# print(av[3])

# a = None
# b = 1
# print(a is not None)
# print(a is None)
# print(b is not None)
