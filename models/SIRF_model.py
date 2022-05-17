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


phase_names = ["0th","1st","2nd","3rd","4th","5th"]


# SIRF subclass
class SIRF_model(SIR_model):
    
    
    def check_params(self,params):
        if params["rho"] == None and params["sigma"] == None:
            warnings.warn("No value for rho and sigma given. Estimation of both with ....")
            # self.fit_rho_and_sigma_function()
        else:
            if params["rho"] == None:
                warnings.warn("No value for rho given. Estimation of rho with ....")
                # self.fit_rho_function()
            elif params["sigma"] == None:
                warnings.warn("No value for sigma given. Estimation of sigma with ....")
                # self.fit_sigma_function()
            
        # make sure all required params are given
        if params["kappa"] == None or params["sigma"] == None:
            raise Warning("kappa and sigma needed for SIR-F model")
        
        
        
    def create_sir(self,params,plot=False):
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
        self.create = True

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

    def create_scenario(self, name, scenario_end_list,rho_constant_list=None,sigma_constant_list=None,kappa_constant_list=None,theta_constant_list=None,plot=False):
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
            self.snl.clear(name=name)
            self.scenario_names.append(name)
            
            #todo: only add main scenario once, so that new scenarios can be added without probelms
            for i,(phase_date, rho_constant, sigma_constant, kappa_constant, theta_constant) in enumerate(zip(scenario_end_list,rho_constant_list,sigma_constant_list,kappa_constant_list,theta_constant_list)): 
                # todo: cant main be renamed to something more meaingful?
                phase_date = pd.to_datetime(phase_date)
                # todo: cant main be renamed to something more meaingful?
                # only add phase if it goes beyond last added phase
                if phase_date > self.last_phase_added:
                    self.snl.add(end_date=phase_date, name="Main")
                    self.last_phase_added = phase_date   
                    
                # adjust original param values
                # always update rho and sigma. If no change is desired than same constant has to be given as input                rho_new = self.snl.get("rho", phase=f"{i}th") * rho_constant
                rho_new = self.snl.get("rho", phase=phase_names[i]) * rho_constant
                sigma_new = self.snl.get("sigma", phase=phase_names[i]) * sigma_constant
                kappa_new = self.snl.get("kappa", phase=phase_names[i]) * kappa_constant
                theta_new = self.snl.get("theta", phase=phase_names[i]) * theta_constant
                
                # Add th i-th phase with the newly calculated params
                self.snl.add(end_date=phase_date, name=name, rho=rho_new,sigma=sigma_new,kappa=kappa_new,theta=theta_new)
            
            # print summary
            print(f"{self.snl.summary()}")
            # get dataframe with Infected, Main, and scenario
            self.phases_plot = self.snl.history(target="Infected",show_figure=False)

            mask = np.array([(pd.to_datetime(self.actual_df.index) >= self.phases_plot.index[0]) & (pd.to_datetime(self.actual_df.index) <= self.phases_plot.index[-1])]).reshape(-1)
            fig, ax = plt.subplots()
            self.phases_plot["Actual"] = self.actual_df.loc[mask]["Infected"].values
            ax.plot(self.phases_plot.index,self.phases_plot["Actual"],label="Acutal")
            ax.plot(self.phases_plot.index,self.phases_plot["Main"],label="Main")
            if plot:
                #todo: plot for all values changed
                #questoin: plotted das automatisch? sollen rt und andere immer geplotted werden?
                _ = self.snl.history(target="Rt")
                # _ = self.snl.history(target="Confirmed")
                _ = self.snl.history(target="rho").head()
                _ = self.snl.history(target="kappa").head()
                _ = self.snl.history(target="theta").head()
            return self.phases_plot,self.snl.describe()