# -*- coding: utf-8 -*-
"""
Created on Thu May 12 19:50:49 2022

@author: User
"""


# imports 
import pandas as pd
import covsirphy as cs
import matplotlib.pyplot as plt
from pprint import pprint
import warnings

#idea: implement __iter__, __str__

class SIR_model():
    
    def __init__(self, country, start_date, end_date):
        """
        Accepted date format: yyyy-mm-dd (e.g. 2020-12-24)
        """
        
        self.country = country
        self.start_date = pd.to_datetime(start_date)
        self.end_date = pd.to_datetime(end_date)
        # load data
        self.load_data()
        # data preparation and cleaning
        self.prep_data()

        self.start_cond()
        
        self.main=False
        
        #return self.df_timerange
                

    def load_data(self):
        """
        Loads up-to-date data from kaggle and stores it as object variable.

        Returns
        -------
        None.
        """
        # data will be downloaded if existing files are older than 24h
        data_loader = cs.DataLoader(directory="kaggle/input")
        # The number of cases and population values in jhu format
        self.jhu_data = data_loader.jhu()
        #cases and deaths whole dataset
        self.total_df = self.jhu_data.total()
        #calculate start conditions given start date and cleaned data
        # days simulated
        self.total_days = (self.end_date - self.start_date).days
        
    def prep_data(self):
        """
        extract clean data and subset for defined country
        Returns
        -------
        None.

        """
        # get rid of rows with missing data
        whole_df = self.jhu_data.cleaned()
        whole_country_df = whole_df.loc[(whole_df["Country"] == self.country)].copy()
        # change time format
        whole_country_df['Date'] = pd.to_datetime(whole_country_df['Date'])
        # store whole cleaned data
        self.whole_country_df = whole_country_df
        
        # define a time range
        mask = (self.whole_country_df['Date'] >= self.start_date) & (self.whole_country_df['Date'] <= self.end_date)
        df_timerange = whole_country_df.loc[mask]
        
        # drop provinces
        province_mask = df_timerange["Province"] == "-"
        df_timerange = df_timerange.loc[province_mask]
        df_timerange.index = pd.RangeIndex(len(df_timerange))
        #store subsetted data
        self.df_timerange = df_timerange
    
    def start_cond(self):
        self.start_fatal_recovered = self.df_timerange["Fatal"][0] + self.df_timerange["Recovered"][0]
        self.start_infected = self.df_timerange["Infected"][0]
        self.start_pop = self.df_timerange["Population"][0]
        # as population size is constant, fatal cases have to be subtracted as well
        self.start_susceptible = self.start_pop - self.start_infected - self.start_fatal_recovered
        self.start_dict = {'Fatal or Recovered': self.start_fatal_recovered, 'Infected': self.start_infected,
                           'Susceptible': self.start_susceptible}
        
        
    def create_sir(self,params={'rho': None, 'sigma': None},plot=False):
        # todo: implement SIR-F model as its own class, not just as parameter (inher) (theta, kappa)
        """
        Creates SIR model and returns resulting data frame

        Parameters
        ----------
        params : TYPE, optional
            DESCRIPTION. The default is {'theta': 0.005, 'kappa': 0.005, 'rho': None, 'sigma': None}.
        sir_f : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        # todo: implement sir f
        # todo: look what params are given and update them
        
        self.check_params(params)
        
        # set model parameters
        self.model = cs.SIR
        self.model.EXAMPLE["param_dict"] = params
        self.model.EXAMPLE['y0_dict'] = self.start_dict
        self.model.EXAMPLE['population'] = self.start_pop
        self.model.EXAMPLE['step_n'] = self.total_days
        
        #todo: do check ich nix meh was passiert, stimmt so?
        # Set tau value and start date of records
        self.example_data = cs.ExampleData(tau=1440, start_date=self.start_date)
        
        # print Model name and parameters
        print(f"created {self.model.NAME}-model with:\nparams:\n{self.model.EXAMPLE['param_dict']}\nstarting conditions:\n\t{self.model.EXAMPLE['y0_dict']}\nsimulating: {self.model.EXAMPLE['step_n']} days")
        

        self.area = {"country": self.country}
        # Add records with SIR model
        self.example_data.add(self.model, **self.area)
        
        # Records with model variables
        df = self.example_data.cleaned()
        
        # prepare the result dataframe
        # todo: make this drop more general. keep what we need not drop to be more flexible
        res_df = df.drop(["Country", "Province", "Confirmed", "ISO3", "Population"], axis = 1)
        
        # add results from the model
        res_df["actual Infected"] = self.df_timerange["Infected"]
        self.res_df = res_df
        
        if plot:
            cs.line_plot(self.res_df.set_index("Date"), title=f"Plot of {self.model.NAME} model", y_integer=True)
        return self.res_df
    
    def check_params(self,params):
        if params["rho"] == None and params["sigma"] == None:
            warnings.warn("No value for rho and sigma given. Estimation of both with ....")
            # self.fit_rho_and_sigma_function()
        elif params["rho"] == None:
            warnings.warn("No value for rho given. Estimation of rho with ....")
            # self.fit_rho_function()
        elif params["sigma"] == None:
            warnings.warn("No value for sigma given. Estimation of sigma with ....")
            # self.fit_sigma_function()
    
    def create_main(self):
        self.main = True
        #todo: only create new snl if no previous snl
        self.snl = cs.Scenario(tau = 1440, **self.area)
        self.snl.register(self.example_data)

        # get the records of the scenario instance
        
        #todo: ist das nötig?
        record_df = self.snl.records(show_plot=False)
        
        # Set 0th phase fro eg from 01Sep2020 to 01Dec2020 with preset parameter values
        self.snl.clear(include_past=True)
        #todo: falls sirf, sicherstellen alle parameter vorhande
        #todo: sicherstellen das alle create_sir vorher aufgerufen wurde!!
        #todo: nicht das gleiche param dict wir params von create_sir?
        
        # past phase
        self.snl.add(end_date=self.end_date, model=self.model, **self.model.EXAMPLE["param_dict"])
        
        
    
    def create_scenario(self, name, scenario_end_list,rho_constant_list=None,sigma_constant_list=None,plot=False):
        """
        adding measures and therefore possible changes to parameters

        Parameters
        ----------
         : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # todo: check that all inputs are lists
        #todo: check that all lists have same length
        #todo: create main has to be called before create scenario
        if not self.main:
            raise Warning("create main has to be called before create scenario")
        
        # Add main scenario
        
        #todo: only add main scenario once, so that new scenarios can be added without probelms
        for i,(phase_date, rho_constant, sigma_constant) in enumerate(zip(scenario_end_list,rho_constant_list,sigma_constant_list)): 
            # todo: cant main be renamed to something more meaingful?
            self.snl.add(end_date=phase_date, name="Main")
        
            # add scenario "name"
            self.snl.clear(name=name)
            # adjust original rho value
            # todo: only adjust those params that are given to this funcion
            if rho_constant != None:
                rho_new = self.snl.get("rho", phase=f"{i}th") * rho_constant
            else:
                rho_new = self.snl.get("rho", phase=f"{i}th")
            if sigma_constant != None:
                sigma_new = self.snl.get("sigma", phase=f"{i}th") * sigma_constant
            else:
                sigma_new = self.snl.get("sigma", phase=f"{i}th")

            # Add th i-th phase with the newly calculated params
            self.snl.add(end_date=phase_date, name=name, rho=rho_new,sigma=sigma_new)
        
        # print summary
        print(f"{self.snl.summary()}")
        if plot:
            #todo: plot for all values changed
            #questoin: plotted das automatisch? sollen rt und andere immer geplotted werden?
            _ = self.snl.history(target="Rt")
            # todo: infected plot in this case is the df, where actual values are missing.
            # .. we either have to insert missing actual values beforehand for that time or impede function from plotting and plot it ourselves
            self.z = infected_plot = self.snl.history(target="Infected",show_plot=False)
            # try to get better plot by creating own plot
            fig, ax = plt.subplots()
            print(infected_plot.columns)
            print(infected_plot)
            ax.plot(infected_plot["Name"],infected_plot["Actual"])
            ax.plot(infected_plot["Name"],infected_plot["Lockdown"])
            ax.plot(infected_plot["Name"],infected_plot["Main"])
            _ = self.snl.history(target="rho").head()
            
    def simluate_scenario(self, name):
        #todo: raise error wenn name kein scenario ist
        _ = self.snl.simulate(name = name)

            
    def get_res_df(self):
        """
        Getter-function for resulting datafram

        Returns self.res_df
        -------
        None.
        """
        try:
            return self.res_df
        #todo: insert correct error
        except:
            raise Warning("No res_df. Run create_sir first")


    def trend_analysis(self):
        """
        Analyses the phases of the wave of infections and estimates their respective parameters.

        Returns
        -------
        Dataframe
        """

        # scenario generation
        srt = cs.Scenario(country=self.country, province=None)
        srt.register(self.jhu_data)

        # time period definition
        srt.timepoints(first_date=self.start_date, last_date=self.end_date, today=self.end_date)

        # phase detection
        _ = srt.trend(show_figure=False)

        # Parameter estimation of the defined model
        # Default value of timeout is 180 sec
        srt.estimate(cs.SIR)

        return srt.summary()


    def get_plot(self):
        #idea: funktion um plots zu returnen
        pass

    
    
def main():
    start_date = '2020-09-01'
    end_date = '2020-12-01'
    country = "Switzerland"
    a = SIR_model(country,start_date,end_date)
    a.create_sir(params={'rho': 0.1, 'sigma': 0.1})
    a.create_main()
    a.create_scenario(name="Lockdown",scenario_end_list=["31Mar2021"],rho_constant_list=[0.5],sigma_constant_list=[1],plot=True)

if __name__ == "__main__":
    main()
    