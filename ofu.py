# extract all multiplets from a given dataset
import numpy as np
from ofu_tools import ang_distance,doublet_average

# Definition of the doublets and multiplets: arXiv:1612.06028v2

class multiplets():
    # data: , delta_t: cut on time difference, delta_ang: cut on angular separation, telescope_term: FoV term in doublet TS, floor: clip values of ang Errors
    def __init__(self,data,delta_t=100.,delta_ang=3.5,telescope_term=0.9,floor=0.2):
        # assuming data that is ordered
        assert (np.sum(np.diff(data['time'])<0.)>0 ), 'data has to be ordered in time!'
        self.data = data
        self._n_multiplets = 0
        self._multiplets = []
        self._delta_t = delta_t/86400. #in MJD
        self._delta_ang = np.radians(delta_ang) #in radians
        self._theta_a = telescope_term
        if floor is not None:
            self._apply_floor(floor)
        # scan for multiplets
        self.__process()


    def _apply_floor(self,floor):
        self.data['angErr'] = np.clip(self.data['angErr'],a_min=floor,a_max=None)
        self._floor = floor

    @property
    def floor(self):
        return self._floor
    @property
    def n_multiplets(self):
        return self._n_multiplets
    @property
    def get_doublets(self):
        return self._multiplets[0]

    def __process(self):
    #find all multiplets with multiplicity m
        last_ntet = False
        multiplicity = 2
        while not last_ntet:
            # find current multiplicity in data?
            time_diff = np.abs(np.diff(self.data['time'], n=multiplicity-1 ) )

            if multiplicity==2: #case of the doublets
                space_diff = ang_distance(self.data[:-(multiplicity-1)] , self.data[(multiplicity-1):])
                total_cut = np.logical_and(np.abs(space_diff<self._delta_ang),time_diff<self._delta_t)
                # defines the two events for each doublet:
                event1 = self.data[:-1][total_cut]
                event2 = self.data[1:][total_cut]
                # can also have doublets that are not subsequent events:
                n_lag = 2
                additional_time_diff = np.abs(np.diff(self.data['time'], n=n_lag ))
                while np.sum( additional_time_diff<self._delta_t ):
                    additional_space_diff = ang_distance(self.data[:-n_lag] , self.data[n_lag:])
                    additional_cut = np.logical_and( additional_time_diff<self._delta_t , additional_space_diff<self._delta_ang )
                    if np.sum(additional_cut)>0:
                        #print('Found lagged doublets with n_lag = ',n_lag)
                        # insert events sorted in time:
                        e1_insert = np.searchsorted(event1['time'],self.data[:-n_lag][additional_cut]['time'])
                        e2_insert = np.searchsorted(event2['time'],self.data[n_lag:][additional_cut]['time'])
                        np.insert(event1,e1_insert,self.data[:-n_lag][additional_cut])
                        np.insert(event1,e2_insert,self.data[n_lag:][additional_cut])
                    n_lag += 1
                    additional_time_diff = np.abs(np.diff(self.data['time'], n=n_lag ))

                # at least one pair fulfills criterion
                if np.sum( np.logical_and( np.abs(time_diff<self._delta_t) , space_diff<self._delta_ang ) )>0:
                    self._n_multiplets = multiplicity
                    doublets=doublet_average(event1, event2 ,self._theta_a)
                    self._multiplets.append(doublets)
                    multiplicity+=1
                    # print('Doublets scanned.')

            else:
                # test all sub-pairs, must all fit the spacial distance criterion
                m = multiplicity-1
                space_diff = ang_distance(self.data[:-m] , self.data[m:])
                space_cut = np.abs(space_diff<self._delta_ang)
                for i in range(m,0,-1):
                    for j in range(i):
                        #print('Testing indices: ',m-i,-i,' ',m-j,-j)
                        if j!=0:
                            space_diff_2 = ang_distance(self.data[m-i:-i] , self.data[m-j:-j])
                        else: # slice [a:-0] would give [] ..
                            space_diff_2 = ang_distance(self.data[m-i:-i] ,self.data[m-j:])
                        space_cut = np.logical_and(space_cut,np.abs(space_diff_2<self._delta_ang) )

                if np.sum( np.logical_and( np.abs(time_diff<self._delta_t) , space_cut ) )>0:
                # at least one pair fulfills criterion
                    self._n_multiplets = multiplicity

                # hiher multiplicity doublets have no TS defined, so far nothing done to them
                    print('Found a (at least one) multiplet with multiplicity {} in the dataset!'.format(multiplicity))
                    multiplicity+=1

                else:
                    last_ntet = True
                    print('Multiplets processed')


    def info(self):
        print(self.data.dtype)
        print(np.min(self.data['angErr']))




    #def get_doublets():

    #def get_multiplets():
