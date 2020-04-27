import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as mm
from scipy.stats import norm


class ResStatInfo:

    # used for tracking class's instances and counts
    Instances = {}
    InstantCount = 0
    pd.set_option("display.max_rows", 10, "display.max_columns", 10)

    @classmethod
    def InstanceCount(cls):
        return len(cls.Instances)

    @classmethod
    def GetInstances(cls):
        return (key for key in cls.Instances.keys())

    @classmethod
    def Get(cls, target):
        return cls.Instances.get(target, None)

    def __init__(self, resid):
        self.resid = resid
        self.resname = None     # will implement later
        self.coord_dset = pd.Series(name='coordination data')
        self.dwell_time_dset = pd.Series(name='dwell time data')
        self.notes = []
        self.Instances[self.resid] = self

    def __repr__(self):
        return "ResStatInfo {} with ".format(self.get_resid())

    def __str__(self):
        return 'ResStatInfo-' + self.get_resid()

    def __lt__(self, other):
        if type(self) != type(other):
            raise Exception('Incompatible argument to __lt__:' + str(other))
        return self.get_resid() < other.get_resid()

    def get_resid(self):
        return self.resid

    def get_resname(self):
        return self.resname

    def get_coord_dset(self):
        return self.coord_dset

    def get_dwell_time_dset(self):
        return self.dwell_time_dset

    def update_res_data(self, sim_info, coord, dwell_time):
        self.coord_dset = self.coord_dset.append(pd.Series([coord], index=[sim_info]))
        self.dwell_time_dset = self.dwell_time_dset.append(pd.Series([dwell_time], index=[sim_info]))

    def add_coord(self, coord, sim_info):
        self.coord_dset = self.coord_dset.append(pd.Series([coord], index=[sim_info]))

    def remove_coord(self, sim_info):
        self.coord_dset = self.coord_dset.drop([sim_info])

    def add_dwell_time(self, dwell_time, sim_info):
        self.dwell_time_dset = self.dwell_time_dset.append(pd.Series([dwell_time], index=[sim_info]))

    def remove_dwell_time(self, dwell_time):
        self.dwell_time_dset = self.dwell_time_dset.drop([dwell_time])

    def combine_coords(self):
        wt_coords = []
        mut_coords = []
        for index in self.coord_dset.keys():
            if index.startswith('wt'):
                wt_coords += self.coord_dset[index]
            else:
                mut_coords += self.coord_dset[index]
        return wt_coords, mut_coords

    def get_pval(self, darray1, darray2, null=0, ttype='two_sided'):
        pe1 = np.mean(darray1)
        pe2 = np.mean(darray2)
        se = mm.sqrt(pe1*(1-pe1)/len(darray1) + pe2*(1-pe2)/len(darray2))
        tstat = abs((pe1-pe2-null)/se)
        if ttype == 'one_sided':
            pval = 1-norm.cdf(tstat)
        else:
            pval = 2*(1-norm.cdf(tstat))
        return pval

    def get_coord_stats(self, key='by_type'):
        stat_df = pd.DataFrame(data=None)
        if key == 'by_type':
            wt, mut = self.combine_coords()
            group_stat = {'wt_group': [len(wt), np.mean(wt), self.is_coord_normal(wt), np.nan],
                          'mut_group': [len(mut), np.mean(mut), self.is_coord_normal(mut), np.nan],
                          'difference': [np.nan, np.mean(wt)-np.mean(mut),
                                         self.is_coord_normal(wt) and self.is_coord_normal(mut),
                                         self.get_pval(wt, mut, ttype='one_sided')]}
            stat_df = pd.DataFrame(data=group_stat, index=['count', 'proportion', 'normality', 'P-value'])
        elif key == 'all':
            coord_stats = {key: [len(value), np.mean(value)] for key, value in self.coord_dset.items()}
            stat_df = pd.DataFrame(data=coord_stats, index=['count', 'proportion'])
        return stat_df

    def get_dwell_time_stats(self):
        index = self.dwell_time_dset.keys()[1]
        return self.dwell_time_dset[index]

    def is_coord_normal(self, darray):
        pe = np.mean(darray)    # point estimate
        size = len(darray)      # sample size
        return size*pe >= 10 and size*(1-pe) >= 10

    def get_note(self):
        for note in self.notes:
            yield note
        # return (note for note in self.notes)

    def add_note(self, note):
        pass


class Parser:

    def __init__(self, filename):
        self.filename = filename
        self.src = None
        self.line = None

    def parse(self):
        sim_info = self.filename.rpartition('_')[0].rpartition('/')[2]
        resid_lst = []
        coord_lst = []
        dwell_time_lst = []
        with open(self.filename) as self.src:
            self.read_next_line()
            while self.line:
                resid_lst.append(self.extract_resid())
                coord_lst.append(self.extract_coord())
                dwell_time_lst.append(self.extract_dwell_time())
                self.read_next_line()
            return resid_lst, sim_info, coord_lst, dwell_time_lst

    def read_next_line(self):
        self.line = self.src.readline().rstrip('\n')

    def extract_resid(self):
        return self.line.split(' ')[0]

    def extract_coord(self):
        # print(len(self.line.split(' ')[501:]))
        return tuple(int(elt) for elt in self.line.split(' ')[501:3000])        # cut of first 500 of equilibration time

    def extract_dwell_time(self):
        blocks = self.line.partition(' ')[2][501:3000].replace(' ', '').split('0')  # cut of first 500 of equilibration time
        return tuple(len(block) for block in blocks if block != '')


def create_res_objects_dict(directory, filenames):
    with open(directory + filenames) as file:
        res_objs = dict()
        filename = file.readline().rstrip('\n')
        while filename:
            resids, sim_info, coords, dwell_times = Parser(rawdatadir + filename).parse()
            for n, resid in enumerate(resids):
                res_objs.setdefault(resid, ResStatInfo(resid))
                res_objs[resid].update_res_data(sim_info, coords[n], dwell_times[n])
            filename = file.readline().rstrip('\n')
        return res_objs


def write_coord_stats(res_objs_dict):
    with open(resultdir + 'coord_stats.csv', 'w') as file:
        file.write(','.join(['residue', 'wt_group [count]', 'wt_group [proportion]', 'wt_group [normal]',
                              'mut_group [count]', 'mut_group [proportion]', 'mu_group [normal]',
                              'diff [proportion]', 'diff [normal]', 'diff [p-value]']) + '\n')
        for resid, res_obj in res_objs_dict.items():
            res_df = res_obj.get_coord_stats(key='by_type')
            file.write(','.join([resid, str(res_df.at['count', 'wt_group']), str(res_df.at['proportion', 'wt_group']),
                                 str(res_df.at['normality', 'wt_group']), str(res_df.at['count', 'mut_group']),
                                 str(res_df.at['proportion', 'mut_group']), str(res_df.at['normality', 'mut_group']),
                                 str(res_df.at['proportion', 'difference']), str(res_df.at['normality', 'difference']),
                                 str(res_df.at['P-value', 'difference'])]) + '\n')


def main():
    # create ResStatInfo objects for all amino acid residues
    residues = create_res_objects_dict(rawdatadir, 'filenames.txt')

    # write statistical analysis results based of receptor type to <coord_stats.csv> file
    write_coord_stats(residues)

    # dwell time analysis
    # transformed = [mm.log2(elt) for elt in residues['431'].get_dwell_time_stats()]
    # plt.hist(residues['431'].get_dwell_time_stats())
    # plt.hist(transformed)
    # plt.show()


if __name__ == '__main__':
    rawdatadir = '../rawdata/'
    resultdir = '../results/'
    main()
