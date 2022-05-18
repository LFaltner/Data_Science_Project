# -*- coding: utf-8 -*-
"""
Created on Sun May 15 21:16:54 2022

@author: User
"""

# imports 
import pandas as pd
import covsirphy as cs
import matplotlib.pyplot as plt
from SIR_model import SIR_model
import numpy as np
from pprint import pprint
import warnings

colors = ["r","g","b","m","black","brown"]


# SIRF subclass
class SIRF_model(SIR_model):
    
    
    def start_cond(self):
        """
        Extract starting conditions from data frame.

        """
        self.start_recovered = self.df_timerange["Recovered"][0]
        self.start_fatal = self.df_timerange["Fatal"][0]
        self.start_infected = self.df_timerange["Infected"][0]
        self.start_pop = self.df_timerange["Population"][0]
        # as population size is constant, fatal cases have to be subtracted as well
        self.start_susceptible = self.start_pop - self.start_infected - self.start_fatal - self.start_recovered
        self.start_dict = {'Recovered': self.start_recovered, 'Infected': self.start_infected,
                           'Susceptible': self.start_susceptible, "Fatal":self.start_fatal}
    
    
    def create_sirf(self,params={},plot=False):
        """
        Creates SIRF model with defined parameters.

        """
        self.create = True


        params = self.check_params(params)

        # set model parameters
        self.model = cs.SIRF
        self.model.EXAMPLE["param_dict"] = params
        self.model.EXAMPLE['y0_dict'] = self.start_dict
        self.model.EXAMPLE['population'] = self.start_pop
        self.model.EXAMPLE['step_n'] = self.total_days

        # Set tau value and start date of records
        self.example_data = cs.ExampleData(tau=1440, start_date=self.start_date)

        # print Model name and parameters
        print(f"created {self.model.NAME}-model with:\nparams:\n\t{self.model.EXAMPLE['param_dict']}\nstarting conditions:\n\t{self.model.EXAMPLE['y0_dict']}\nsimulating: {self.model.EXAMPLE['step_n']} days")

        self.area = {"country": self.country}
        # Add records with SIR model
        self.example_data.add(self.model, **self.area)

        # Records with model variables
        df = self.example_data.cleaned()

        # prepare the result dataframe
        res_df = df.drop(["Country", "Province", "Confirmed", "ISO3", "Population"], axis=1)

        # add results from the model
        res_df["actual Infected"] = self.df_timerange["Infected"]
        self.res_df = res_df

        if plot:
            cs.line_plot(self.res_df.set_index("Date"), title=f"Plot of {self.model.NAME} model", y_integer=True)
        return self.res_df
    
    def check_params(self,params):
        """
        Check for correctness of input params.
        """
        
        # make sure all required params are given
        if "kappa" not in params.keys() or "theta" not in params.keys():
            raise Warning("kappa and sigma needed for SIR-F model")
            
        if "rho" not in params.keys() and "simga" not in params.keys():
            print("No values for rho and sigma given. Estimation of both parameters, this may take a while...")
            params["rho"], params["sigma"] = self.parameter_estimation()
            print(f"Estimations: Rho = {params['rho']}, Sigma = {params['sigma']}\n\n")
            return params
        
        elif "rho" not in params.keys():
            print("No value for rho given. Estimation of rho, this may take a while...")
            params["rho"], _ = self.parameter_estimation()
            print(f"Estimation: Rho = {params['rho']}\n\n")
            return params
        
        elif "sigma" not in params.keys():
            print("No value for sigma given. Estimation of sigma, this may take a while...")
            _, params["sigma"] = self.parameter_estimation()
            print(f"Estimation: Sigma = {params['sigma']}\n\n")
            return params
        else:
            return params
            
        

    def create_scenario(self, name, scenario_end_list,rho_constant_list=None,sigma_constant_list=None,kappa_constant_list=None,theta_constant_list=None,plot=False):
        """
        Creating scenarios: gives ability to change parameters at specific points in times, which faciliates the modeling of real world measures.

        """

        if not self.main:
            raise Warning("create main has to be called before create scenario")
        if name not in self.scenario_names:
            self.scenario_names.append(name)
            # self.colors[name] = self.colors_list[len(self.scenario_names)-1]
        else:
            raise Warning("Scenario already exists")
        
        # Add main scenario
        self.snl.clear(name=name)
        

        for i,(phase_date, rho_constant, sigma_constant, kappa_constant, theta_constant) in enumerate(zip(scenario_end_list,rho_constant_list,sigma_constant_list,kappa_constant_list,theta_constant_list)): 
            phase_date = pd.to_datetime(phase_date)
            # only add phase if it goes beyond last added phase
            if phase_date > self.last_phase_added:
                self.snl.add(end_date=phase_date, name="Main")
                self.last_phase_added = phase_date   
                
            # adjust original param values
            # always update rho and sigma. If no change is desired than same constant has to be given as input                rho_new = self.snl.get("rho", phase=f"{i}th") * rho_constant
            rho_new = self.snl.get("rho", phase="0th") * rho_constant
            sigma_new = self.snl.get("sigma", phase="0th") * sigma_constant
            kappa_new = self.snl.get("kappa", phase="0th") * kappa_constant
            theta_new = self.snl.get("theta", phase="0th") * theta_constant
            
            # Add th i-th phase with the newly calculated params
            self.snl.add(end_date=phase_date, name=name, rho=rho_new,sigma=sigma_new,kappa=kappa_new,theta=theta_new)
        
        # print summary
        print(f"{self.snl.summary()}")
        # get dataframe with Infected, Main, and scenario
        self.infected_plot = self.snl.history(target="Infected",show_figure=False)
        # get dataframe for total confirmed cases
        self.confirmed_plot = self.snl.history(target="Confirmed",show_figure=False)

        mask = np.array([(pd.to_datetime(self.actual_df.index) >= self.infected_plot.index[0]) & (pd.to_datetime(self.actual_df.index) <= self.infected_plot.index[-1])]).reshape(-1)
        self.infected_plot["Actual"] = self.actual_df.loc[mask]["Infected"].values
        self.confirmed_plot["Actual"] = self.actual_df.loc[mask]["Confirmed"].values
        
        # prepare plots
        fig, (ax_inf,ax_conf) = plt.subplots(nrows=2,ncols=1,sharex=True)
        ax_inf.set_title("Active Infected")
        ax_conf.set_title("Total confirmed cases")
        # plot actual and main infected
        ax_inf.plot(self.infected_plot.index,self.infected_plot["Actual"],label="Acutal")
        ax_inf.plot(self.infected_plot.index,self.infected_plot["Main"],label="Main")
        # plot actual and main confirmed
        ax_conf.plot(self.confirmed_plot.index,self.confirmed_plot["Actual"],label="Acutal")
        ax_conf.plot(self.confirmed_plot.index,self.confirmed_plot["Main"],label="Main")
        
        
        for name in self.scenario_names:
            ax_inf.plot(self.infected_plot.index,self.infected_plot[name],label=name)
            ax_conf.plot(self.confirmed_plot.index,self.confirmed_plot[name],label=name)
        plt.legend()
        
        if plot:
            _ = self.snl.history(target="Rt")
            _ = self.snl.history(target="rho")
            _ = self.snl.history(target="sigma")
            _ = self.snl.history(target="kappa")
            _ = self.snl.history(target="theta")

                