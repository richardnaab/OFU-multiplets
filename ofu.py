# extract all multiplets from a given dataset
import numpy as np
from ofu_tools import ofu_ang_distance,ofu_doublet_average


def sub_diff(arr,n):
    return arr[n:]-arr[:-n]

def n_angular_pairs(evts):
    n_mult = 0
    for count in range(1,len(evts)):
        space_diff = np.abs( ofu_ang_distance(evts[:-count] , evts[count:]) )
        n_mult += np.sum(space_diff<self._delta_ang)
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
            print('Found any doublet pair?',np.any(total_cut))
            if np.any(total_cut):
                # defines the two events for each doublet:
                if first_events:
                    event1 = self.data[:-n_lag][total_cut]
                    event2 = self.data[n_lag:][total_cut]
                    first_events = False
                else: # insert at correct times
                    e1_insert = np.searchsorted(event1['time'],self.data[:-n_lag][total_cut]['time'])
                    e2_insert = np.searchsorted(event2['time'],self.data[n_lag:][total_cut]['time'])
                    np.insert(event1,e1_insert,self.data[:-n_lag][total_cut])
                    np.insert(event1,e2_insert,self.data[n_lag:][total_cut])
            n_lag += 1
            time_diff = sub_diff(self.data['time'],n_lag)

        # at least one pair fulfills criterion
        if len(event1)>0:
            doublets=ofu_doublet_average(event1, event2 ,self._theta_a)
            self._n_multiplets = 2
            self._multiplets.append(doublets)
            #print('Doublets scanned.')

            # scan here for higher multiplicities:
            uni,index,inv,counts = np.unique(event2,return_index=True, return_inverse=True, return_counts=True)
            if np.any(counts>=2):
                #print('Test for higher multiplicities:')
                candidates, = np.where(counts>=2)
                for cand in candidates:
                    mask = inv==cand    # event1[mask] are all events that form a doublet with cand
                    # find number of pair wise doublets (so far only to check if any triplet!)
                    # -> still head for improvement, e.g. really save triplets or higher multiplicity multiplets
                    if n_angular_pairs(event1[mask]) > 0:
                        self._n_multiplets = 3
                        print('Found {} tiplets in the sample!'.format(n_angular_pairs(event1[mask])) )

        else:
            print('No doublets found')


        # # can also have doublets that are not subsequent events:
        # n_lag = 2
        # additional_time_diff = np.abs(np.diff(self.data['time'], n=n_lag ))
        # while np.sum( additional_time_diff<self._delta_t ):
        #     additional_space_diff = ofu_ang_distance(self.data[:-n_lag] , self.data[n_lag:])
        #     additional_cut = np.logical_and( additional_time_diff<self._delta_t , additional_space_diff<self._delta_ang )
        #     if np.sum(additional_cut)>0:
        #         #print('Found lagged doublets with n_lag = ',n_lag)
        #         # insert return self._multiplets[0]events sorted in time:
        #         e1_insert = np.searchsorted(event1['time'],self.data[:-n_lag][additional_cut]['time'])
        #         e2_insert = np.searchsorted(event2['time'],self.data[n_lag:][additional_cut]['time'])
        #         np.insert(event1,e1_insert,self.data[:-n_lag][additional_cut])
        #         np.insert(event1,e2_insert,self.data[n_lag:][additional_cut])
        #     n_lag += 1
        #     additional_time_diff = np.abs(np.diff(self.data['time'], n=n_lag ))


            # last_ntet = False
            # multiplicity = 2
            # while not last_ntet:
            #     # find current multiplicity in data?
            #     #time_diff = np.abs(np.diff(self.data['time'], n=multiplicity-1 ) )
            #
            #     if multiplicity==2: #case of the doublets

            # else:
            #     if self._only_doublets:
            #         last_ntet = True
            #         print('Only Doublets processed')
            #     else:
            #         #time_diff = np.abs(np.diff(self.data['time'], n=multiplicity-1 ) )
            #         time_diff = sub_diff(self.data['time'],multiplicity-1)
            #         # test all sub-pairs, where at least x pairs have to fulfill the spacial distance criterion
            #         m = multiplicity-1
            #         space_diff = np.abs( ofu_ang_distance(self.data[:-m] , self.data[m:]) )
            #         space_cut = space_diff<self._delta_ang
            #         for i in range(m,0,-1):
            #             for j in range(i):
            #                 #print('Testing indices: ',m-i,-i,' ',m-j,-j)
            #                 if j!=0:
            #                     space_diff_2 = ofu_ang_distance(self.data[m-i:-i] , self.data[m-j:-j])
            #                 else: # slice [a:-0] would give [] ..
            #                     space_diff_2 = ofu_ang_distance(self.data[m-i:-i] ,self.data[m-j:])
            #                 space_cut = np.logical_and(space_cut,np.abs(space_diff_2<self._delta_ang) )
            #
            #         if np.sum( np.logical_and( np.abs(time_diff<self._delta_t) , space_cut ) )>0:
            #         # at least one pair fulfills criterion
            #             self._n_multiplets = multiplicity
            #
            #         # hiher multiplicity doublets have no TS defined, so far nothing done to them
            #             print('Found a (at least one) multiplet with multiplicity {} in the dataset!'.format(multiplicity))
            #             multiplicity+=1
            #
            #         else:
            #             last_ntet = True
            #             print('Multiplets processed')


    def info(self):
        print(self.data.dtype)
        print(np.min(self.data['angErr']))




    #def get_doublets():

    #def get_multiplets():
