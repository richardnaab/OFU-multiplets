# extract all multiplets from a given dataset
import numpy as np
from ofu_tools import ofu_ang_distance,ofu_doublet_average,ofu_multiplet_average, multiplet_dtype,restrict_trig_radius


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
    # data: , delta_t: cut on time difference (sec), delta_ang: cut on angular separation (deg),
    # telescope_term: FoV term (deg) in doublet TS, floor: clip values of ang Errors (deg),
    # info: some prints for help, best_trig: if multiple doublets have same trigger, save only the one with highest TS?
    # triplet_crit: criterion for the triplets to be selected. "loose" requires one event to belong to (at least) two doublets,
    #                                                       "strict" requires the two remaining events to also form a doublet
    def __init__(self,data,delta_t=100.,delta_ang=3.5,telescope_term=0.9,floor=0.2,info=False,best_trig=True,triplet_crit='loose'):
        # assuming data that is ordered
        assert (np.sum(np.diff(data['time'])<0.)>=0 ), 'data has to be ordered in time!'
        self.data = data
        self._n_multiplets = 0
        self._multiplets = []
        self._delta_t = delta_t/86400. #in MJD
        self._delta_ang = np.radians(delta_ang) #in radians
        self._theta_a = telescope_term
        self._info = info
        self._best_trig = best_trig
        if floor is not None:
            self._apply_floor(floor)
        if triplet_crit in {'loose','strict'}:
            self._triplet_crit = triplet_crit
            # scan for multiplets
            self.__process()
        else:
            print('Triplet criterion not understood, nothing done!')

    def _apply_floor(self,floor):
        self.data['angErr'] = np.clip(self.data['angErr'],a_min=np.radians(floor),a_max=None)
        self._floor = floor

    def shuffle_ra(self):
        self.data['ra'] = np.random.permutation(self.data['ra'])

    def triggers_in_radius(self,ra,dec,radius):
        # coordinates and radius in radians
        if self._n_multiplets >=2:
            self._multiplets[0] = restrict_trig_radius(self._multiplets[0],ra,dec,radius)

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
                else:
                    event1 = np.append(event1,self.data[:-n_lag][total_cut])
                    event2 = np.append(event2,self.data[n_lag:][total_cut])
            n_lag += 1
            time_diff = sub_diff(self.data['time'],n_lag)

        # at least one pair fulfills criterion
        if not first_events:
            self._n_multiplets = 2
            doublets=ofu_doublet_average(event1, event2 ,self._theta_a)
            # scan here for higher multiplicities:
            multiplets = np.zeros(0, dtype=multiplet_dtype)
            # find triplets/quadruplets/..
            all_times = np.append(event1['time'],event2['time'])
            uni,index,counts = np.unique(all_times, return_index=True,return_counts=True)

            if self._best_trig :
                # find doublets in which both events belong to one doublet only
                unique_ind = index[counts==1]
                ev1_uni_ind = unique_ind[unique_ind<len(all_times)/2]
                ev2_uni_ind = unique_ind[unique_ind>=len(all_times)/2]-len(all_times)/2
                # ""-len(all_times)/2" shifts possible indices back to [0,....,len(event2)-1]
                ev1_uni_bool = np.full(len(all_times)/2,False)
                ev1_uni_bool[ev1_uni_ind] = True
                ev2_uni_bool = np.full(len(all_times)/2,False)
                ev2_uni_bool[ev2_uni_ind] = True
                #doublets_best_trig = doublets[np.logical_and(ev1_uni_bool,ev2_uni_bool)]
                multiplets = np.append(multiplets,doublets[np.logical_and(ev1_uni_bool,ev2_uni_bool)])

            else:
                multiplets = np.append(multiplets,doublets)

            # Distinguish the triplet criterions:
            if self._triplet_crit == 'loose':
                # go through all events that belong to at least two doublets:
                if np.any(counts>=2):
                    # at least one event belongs to at least two doublets
                    self._n_multiplets = 3
                    if self._info:
                        print('Found {} triplets (loose criterion) in the sample!'.format(np.sum(counts>=2)) )
                for cand in uni[counts>=2]:
                    #select all events contributing to the multiplet:
                    central = self.data[np.where(self.data['time']==cand)]
                    before = event1[np.where(event2['time']==cand)]
                    after = event2[np.where(event1['time']==cand)]

                    mult = np.append(central,before)
                    mult = np.append(mult,after)
                    #calculate combined position: (no TS in this case)
                    multiplets = np.append(multiplets,ofu_multiplet_average(mult,3))

            elif self._triplet_crit == 'strict':
                # all doublets with both events belonging to one doublet only are found.
                # array to store pair-wise doublets that pass the loose triplet criterion but fail the strict one
                sub_doublets = np.zeros(0, dtype=multiplet_dtype)
                higher_mult = np.zeros(0, dtype=multiplet_dtype)

                # go through all events that belong to at least two doublets
                for cand in uni[counts>=2]:
                    #print('Next triplet cand')
                    #select all events contributing to the multiplet:
                    central = self.data[np.where(self.data['time']==cand)]
                    before = event1[np.where(event2['time']==cand)]
                    after = event2[np.where(event1['time']==cand)]

                    # need to check: temporal AND spacial overlap, again!
                    cluster_events = np.append(before,central)
                    cluster_events = np.append(cluster_events,after)
                    cluster_times = np.sort(cluster_events['time'])
                    # Loop over all possible time combinations:
                    i=0
                    test_time = cluster_times[i]
                    evt_in_cluster = cluster_events[np.logical_and(cluster_events['time']>=test_time,cluster_events['time']<=test_time+self._delta_t)]
                    # catch case of triplet .. higher multiplicity in any case gives triplet, but can miss some of the sub_doublets
                    if len(cluster_events)==3:
                        if len(evt_in_cluster)<3:
                            sub_doublets = np.append(sub_doublets,doublets[doublets['t0']==cand])
                            sub_doublets = np.append(sub_doublets,doublets[doublets['t1']==cand])

                    while len(evt_in_cluster)>=3:
                        print(len(evt_in_cluster))
                        # check if events in cluster also form a doublet with each other
                        if n_angular_pairs(evt_in_cluster,self._delta_ang) > 0:
                            self._n_multiplets = 3
                            #trigger is "central" event (corresponding to the cand time)
                            higher_mult = np.append(higher_mult,ofu_multiplet_average(evt_in_cluster,3))
                        else:
                            pass
                            ### To be further completed to also collect sub_doublets, only chose best fit of these ..
                        i += 1
                        test_time = cluster_times[i]
                        evt_in_cluster = cluster_events[np.logical_and(cluster_events['time']>=test_time,cluster_events['time']<=test_time+self._delta_t)]

                higher_mult = np.unique(higher_mult)
                # go over the remaining sub_doublets:
                sub_doublets = np.unique(sub_doublets)
                multiplets = np.append(multiplets,sub_doublets)
                multiplets = np.append(multiplets,higher_mult)


            # at the moment: double counting of triplets ->
            # if triplet fulfills strict criterion, different cand can give the same multiplett
            # keeping them allows to select trigger position ..

            self._multiplets.append(multiplets)

        else:
            if self._info:
                print('No doublets found')


    def info(self):
        print(self.data.dtype)
        print(np.min(self.data['angErr']))


    #def get_multiplets():
