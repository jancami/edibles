import pandas as pd
from edibles.edibles_settings import edibles_pythondir


class FilterDR(object):
    def __init__(self, init_df=None):
        if init_df is None:
            self.df = pd.read_csv(edibles_pythondir + '/data/DR3_ObsLog.csv')
        else:
            self.df = init_df

    def reset_index(func):
        """
        when you put in this decorator, it will automatically return an
        instance of the object and reset the index
        
        allows this:
        FilterDR().filterStar('HD170740').filterOrder(order=[12, 13])
        """

        def reset(self, *args, **kwargs):
            func(self, *args, **kwargs)
            self.df = self.df.reset_index(drop=True)
            return self

        return reset

    def __str__(self):
        pd.set_option('max_colwidth', -1)
        return self.df.to_string()

    @staticmethod
    def parse_time(timestamp):
        """
        :param timestamp: ex. 2014-10-29T07:01:33.557
        :return: 20141029
        """
        return timestamp.split('T')[0].replace('-', '')

    @reset_index
    def filterAll(self, **kwargs):
        """
        FilterDR().filterAll(star=star, date=time, wavelength=lab_wavelength)
            is the same as
        FilterDR().filterStar(star).filterDate(time).filterRange(lab_wavelength).filterOrder()

        INPUT:          type
            star:       [string]                    ex. 'HD170740'
            date:       [string]                    ex. '20150626',
            wavelength: [list of float] or [float]  Filters out files that do not contain any wavelengths in list
            order:      [list of int] or [int]      11 for _O11
            combined:   [boolean]                   If filename contains L, R or U
        """
        keys = kwargs.keys()
        if 'star' in keys:
            self.filterStar(kwargs['star'])
        if 'date' in keys:
            self.filterDate(kwargs['date'])
        if 'wavelength' in keys:
            self.filterRange(kwargs['wavelength'])
        if 'order' in keys:
            self.filterOrder(order=kwargs['order'])
        elif 'combined' in keys:
            self.filterOrder(combined=kwargs['combined'])
        else:
            self.filterOrder()


    @reset_index
    def filterStar(self, star):
        self.df = self.df[self.df.Object == star]

    @reset_index
    def filterRange(self, lab_wavelength):
        if type(lab_wavelength) == list:
            wave1 = lab_wavelength[0]
            filt = (self.df.WaveMax > wave1) & (self.df.WaveMin < wave1)
            for wave in lab_wavelength[1:]:
                filt = filt | ((self.df.WaveMax > wave) & (self.df.WaveMin < wave))
            self.df = self.df[filt]
        else:
            self.df = self.df[(self.df.WaveMax > lab_wavelength) & (self.df.WaveMin < lab_wavelength)]

    @reset_index
    def filterOrder(self, order=[], combined=False):
        assert not (combined and len(order) > 0), "Contradicting inputs"
        if type(order) == int:
            order = [order]
        order = [str(o) for o in order]
        # bit convoluted but it works
        if len(order) > 0:
            self.df = self.df[
                self.df.Filename == self.df.Filename.apply(lambda x: x if x.split('O')[-1][:-5] in order else None)]
        elif combined:
            self.df = self.df[~self.df.Filename.str.contains('O')]
        else:
            self.df = self.df[self.df.Filename.str.contains('O')]

    @reset_index
    def filterDate(self, date):
        df2 = self.df.copy()
        df2.DateObs = df2.DateObs.apply(self.parse_time)  # format date
        self.df = self.df[(df2.DateObs == date)]

    def getCopy(self):
        """
        Returns copy of FilterDR object.
        """
        return FilterDR(self.df)

    def getAllFileNames(self):
        return list(self.df.Filename)

    def getDataFrame(self):
        return self.df

    def getOrders(self):
        """
        Returns unique list of ints of all orders in list

        ex.
        filter = FilterDR()
        filter.filterOrder(order=[12, 13])
        print(filter.getOrders())

        OUTPUT: 
        [12, 13]
        """
        return [int(x) for x in list(set(self.df.Filename.apply(lambda x: x.split('O')[-1][:-5]))) if '/' not in x]

    def getDates(self):
        """
        Returns unique list of strings for all dates (parsed)

        ex.
        filter = FilterDR()
        filter.filterStar('HD170740')
        print(filter.getDates())

        OUTPUT: 
        ['20150626', '20170701', '20150424', '20160613', '20140916', '20170705', '20160505', '20160612', '20140915']
        """
        return list(set(self.df.DateObs.apply(self.parse_time)))

    def getStars(self):
        """
        Returns list of stars

        ex.
        filter = FilterDR()
        print(filter.getStars())

        OUTPUT:
        ['HD116852', 'HD66194', 'HD147889', ... , 'HD 148937', 'HD 133518', 'HD36822']
        """
        return list(set(self.df.Object))


if __name__ == '__main__':
    filter = FilterDR()
    #print(filter.getStars())
    filter.filterStar('HD170740')
    filter.filterOrder(order=[12, 13])
    # filter.filterOrder()
    # filter.filterDate('20140915')
    # filter.filterRange([3300, 5890])
    #print(filter)
    print(filter.getDates())
