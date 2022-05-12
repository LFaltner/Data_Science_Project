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
    
    def __init(self, country, start_date, end_date):
        # todo: define type and other requirments (e.g date format, making sure all needed columns) for variables
        self.country = country
        self.start_date = start_date
        self.end_date = end_date
        # load data
        self.load_data()
        # data preparation and cleaning

                

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
        # todo: make sure date in dataset: raise warning
        self.start_cond()
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
        whole_df = whole_df.loc[(whole_df["Country"] == self.country)].copy()
        # change time format
        whole_df['Date'] = pd.to_datetime(whole_df['Date'])
        # store whole cleaned data
        self.whole_df = whole_df
        
        # define a time range
        # todo: when has date format to be changed
        start_date = pd.to_datetime(self.start_date)
        end_date = pd.to_datetime(self.end_date)
        mask = (whole_df['Date'] >= start_date) & (whole_df['Date'] <= end_date)
        df_timerange = whole_df.loc[mask]
        
        # todo: is this province shit really needed??
        # drop provinces
        province_mask = df_timerange["Province"] == "-"
        df_timerange = df_timerange.loc[province_mask]
        df_timerange.index = pd.RangeIndex(len(df_timerange))
        #store subsetted data
        self.df_timerange = df_timerange
        # todo: return/print df_timerange?
    
    def start_cond(self):
        self.start_fatal_recovered = self.df_timerange["Fatal"][0] + self.df_timerange["Recovered"][0]
        self.start_infected = self.df_timerange["Indecter"][0]
        self.start_pop = self.df_timerange["Population"][0]
        self.start_susceptible = self.start_pop - self.start_infected - self.start_fatal_recovered
        
        self.start_dict = {'Fatal or Recovered': self.start_fatal_recovered, 'Infected': self.start_infected,
                           'Susceptible': self.start_susceptible}
        
        
    def create_sir(self,params={'theta': 0.005, 'kappa': 0.005, 'rho': None, 'sigma': None},sir_f = False,plot=False):
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
        
        if params["rho"] == None and params["sigma"] == None:
            warnings.warn("No value for rho and sigma give. Estimation of both with ....")
            # self.fit_rho_and_sigma_function()
        elif params["rho"] == None:
            warnings.warn("No value for rho give. Estimation of rho with ....")
            # self.fit_rho_function()
        elif params["sigma"] == None:
            warnings.warn("No value for sigma give. Estimation rho with ....")
            # self.fit_sigma_function()
        
        # set model parameters
        model = cs.SIR
        model.EXAMPLE["param_dict"] = params
        model.EXAMPLE['y0_dict'] = self.start_dict
        model.EXAMPLE['population'] = self.pop
        model.EXAMPLE['step_n'] = self.total_days
        
        #todo: do check ich nix meh was passiert, stimmt so?
        # Set tau value and start date of records
        # todo: date format richtig?
        self.example_data = cs.ExampleData(tau=1440, start_date=self.start_date)
        
        # print Model name and parameters
        print(f"created {cs.SIR.NAME}-model with:\nparams: {cs.SIR.EXAMPLE['param_dict']}\nstarting conditions:\n\t{cs.SIR.EXAMPLE['y0_dict']}\nsimulating: {cs.SIR.EXAMPLE['step_n']} days")
        
        # todo: wieso neue instance von cs.SIR???
        model = cs.SIR
        self.area = {"country": self.country}
        # Add records with SIR model
        example_data.add(model, **area)
        
        # Records with model variables
        df = example_data.cleaned()
        
        # prepare the result dataframe
        # todo: make this drop more general. keep what we need not drop to be more flexible
        res_df = df.drop(["Country", "Province", "Confirmed", "ISO3", "Population"], axis = 1)
        
        # add results from the model
        res_df["actual Infected"] = self.df_timerange["Infected"]
        self.res_df = res_df
        
        if plot:
            #todo: how is it possible that more recovered than infected?
            cs.line_plot(self.res_df.set_index("Date"), title=f"Plot of {model.NAME} model", y_integer=True)
        return self.res_df
    
    def create_scenario(self, name, scenario_end,rho_constant=None,sigma_constant=None,kappa_constant=None,theta_constant=None,plot=False):
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
        #todo: only create new snl if no previous snl
        self.snl = cs.Scenario(tau = 1440, **self.area)
        self.snl.register(self.example_data)

        # get the records of the scenario instance
        #todo: ist das nötig?
        record_df = self.snl.records(show_plot=False)
        
        # Set 0th phase fro eg from 01Sep2020 to 01Dec2020 with preset parameter values
        self.snl.clear(include_past=True)
        #todo: time format?
        #todo: hier automatisch sirf? geht auch sir
        #todo: falls sirf, sicherstellen alle parameter vorhande
        #todo: sicherstellen das alle create_sir vorher aufgerufen wurde!!
        #todo: nicht das gleiche param dict wir params von create_sir?
        
        self.snl.add(end_date=self.end_date, model=cs.SIRF, **cs.SIR.EXAMPLE["param_dict"])
        # Add main scenario
        #todo: only add main scenario once, so that new scenarios can be added without probelms
        self.snl.add(end_date=scenario_end, name="Main (normal sir)")
        
        # add scenario "name"
        self.snl.clear(name="Lockdown")
        # adjust original rho value
        # todo: only adjust those params that are given to this funcion
        if rho_constant != None:
            rho_new = self.snl.get("rho", phase="0th") * rho_constant
        else:
            rho_new = self.snl.get("rho", phase="0th")
        if sigma_constant != None:
            sigma_new = self.snl.get("sigma", phase="0th") * sigma_constant
        else:
            sigma_new = self.snl.get("kappa", phase="0th")
        if kappa_constant != None:
            kappa_new = self.snl.get("kappa", phase="0th") * kappa_constant
        else:
            kappa_new = self.snl.get("kappa", phase="0th")
        if theta_constant != None:
            theta_new = self.snl.get("theta", phase="0th") * theta_constant
        else:
            theta_new = self.snl.get("theta", phase="0th")
        # Add th 1st phase with the newly calculated params
        self.snl.add(end_date=scenario_end, name=name, rho=rho_new, kappa=kappa_new, sigma=sigma_new)
        
        # print summary
        print(f"{self.snl.summary()}")
        if plot:
            #todo: plot for all values changed
            #questoin: plotted das automatisch? sollen rt und andere immer geplotted werden?
            _ = self.snl.history(target="Rt")
            _ = self.snl.history(target="Infected")
            self.snl.history(target="rho").head()
            
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

    
    def get_plot(self):
        #idea: funktion um plots zu returnen
        pass
    
    
    
    
    
    
    
    
    
    