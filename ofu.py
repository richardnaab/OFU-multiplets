# extract all multiplets from a given dataset
import numpy as np
from ofu_tools import ofu_ang_distance,ofu_doublet_average


def sub_diff(arr,n):
    return arr[n:]-arr[:-n]

def n_angular_pairs(evts,ang_dist):
    n_mult = 0
    for count in range(1,len(evts)):
        space_diff = np.abs( ofu_ang_distance(evts[:-count] , evts[count:]) )
        n_mult += np.sum(space_diff<ang_dist)
    return n_mult

# Definition of the doublets and multiplets: arXiv:1612.06028v2

class multiplets():
    # data: , delta_t: cut on time difference, delta_ang: cut on angular separation, telescope_term: FoV term in doublet TS, floor: clip values of ang Errors
    def __init__(self,data,delta_t=100.,delta_ang=3.5,telescope_term=0.9,floor=0.2,only_doublets=False):
        # assuming data that is ordered
        assert (np.sum(np.diff(data['time'])<0.)>=0 ), 'data has to be ordered in time!'
        self.data = data
        self._n_multiplets = 0
        self._multiplets = []
        self._delta_t = delta_t/86400. #in MJD
        self._delta_ang = np.radians(delta_ang) #in radians
        self._theta_a = telescope_term
        self._only_doublets = only_doublets
        if floor is not None:
            self._apply_floor(floor)
        # scan for multiplets
        self.__process()


    def _apply_floor(self,floor):
        self.data['angErr'] = np.clip(self.data['angErr'],a_min=np.radians(floor),a_max=None)
        self._floor = floor

    def shuffle_ra(self):
        self.data['ra'] = np.random.permutation(self.data['ra'])

    @property
    def floor(self):
        return self._floor
    @property
    def n_multiplets(self):
        return self._n_multiplets
    @property
    def get_doublets(self):
        if self._n_multiplets >=2:
            return self._multiplets[0]
        else:
            return None

    def __process(self):
    #find all multiplets with multiplicity m
        first_events = True
        # test all pairs, n_lag counting the index difference of tested pair
        n_lag = 1
        time_diff = sub_diff(self.data['time'],n_lag)

        while np.any( time_diff<self._delta_t ):
            space_diff = np.abs( ofu_ang_distance(self.data[:-n_lag] , self.data[n_lag:]) )
            total_cut = np.logical_and(space_diff<self._delta_ang,time_diff<self._delta_t)
            if np.any(total_cut):
                # defines the two events for each doublet: (event 2 is the "trigger" in realtime case)
                if first_events == True:
                    event1 = np.copy(self.data[:-n_lag][total_cut])
                    event2 = np.copy(self.data[n_lag:][total_cut])
                    first_events = False
                else: # insert at correct times
                    #e1_insert = np.searchsorted(event1['time'],self.data[:-n_lag][total_cut]['time'])
                    #e2_insert = np.searchsorted(event2['time'],self.data[n_lag:][total_cut]['time'])
                    #event1 = np.insert(event1,e1_insert,self.data[:-n_lag][total_cut])
                    #event2 = np.insert(event2,e2_insert,self.data[n_lag:][total_cut])
                    event1 = np.append(event1,self.data[:-n_lag][total_cut])
                    event2 = np.append(event2,self.data[n_lag:][total_cut])
            n_lag += 1
            time_diff = sub_diff(self.data['time'],n_lag)

        # at least one pair fulfills criterion
        if not first_events:
            doublets=ofu_doublet_average(event1, event2 ,self._theta_a)
            self._n_multiplets = 2
            self._multiplets.append(doublets)

            # scan here for higher multiplicities:
            uni,index,inv,counts = np.unique(event2,return_index=True, return_inverse=True, return_counts=True)
            if np.any(counts>=2):
                candidates, = np.where(counts>=2)
                for cand in candidates:
                    mask = inv==cand    # event1[mask] are all events that form a doublet with cand
                    # find number of pair wise doublets (so far only to check if any triplet!)
                    # -> still head for improvement, e.g. really save triplets or higher multiplicity multiplets
                    if n_angular_pairs(event1[mask],self._delta_ang) > 0:
                        self._n_multiplets = 3
                        print('Found {} triplets in the sample!'.format(n_angular_pairs(event1[mask],self._delta_ang)) )

        else:
            print('No doublets found')


    def info(self):
        print(self.data.dtype)
        print(np.min(self.data['angErr']))


    #def get_multiplets():
