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
    
    def create_sirf(self,params={},plot=False):
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

        params = self.check_params(params)

        # set model parameters
        self.model = cs.SIRF
        self.model.EXAMPLE["param_dict"] = params
        self.model.EXAMPLE['y0_dict'] = self.start_dict
        self.model.EXAMPLE['population'] = self.start_pop
        self.model.EXAMPLE['step_n'] = self.total_days

        # todo: do check ich nix meh was passiert, stimmt so?
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
        # todo: make this drop more general. keep what we need not drop to be more flexible
        res_df = df.drop(["Country", "Province", "Confirmed", "ISO3", "Population"], axis=1)

        # add results from the model
        res_df["actual Infected"] = self.df_timerange["Infected"]
        self.res_df = res_df

        if plot:
            cs.line_plot(self.res_df.set_index("Date"), title=f"Plot of {self.model.NAME} model", y_integer=True)
        return self.res_df
    
    def check_params(self,params):
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
            if name not in self.scenario_names:
                self.scenario_names.append(name)
            else:
                raise Warning("Scenario already exists")
            
            # Add main scenario
            self.snl.clear(name=name)
            
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
            # todo: make sure smae scenario has same color in all plots
            # todo: plot estimated "rho" of real life?
            # get dataframe with Infected, Main, and scenario
            self.infected_plot = self.snl.history(target="Infected",show_figure=False)
            # get dataframe for total confirmed cases
            self.confirmed_plot = self.snl.history(target="Confirmed",show_figure=False)
    
            mask = np.array([(pd.to_datetime(self.actual_df.index) >= self.infected_plot.index[0]) & (pd.to_datetime(self.actual_df.index) <= self.infected_plot.index[-1])]).reshape(-1)
            self.infected_plot["Actual"] = self.actual_df.loc[mask]["Infected"].values
            self.confirmed_plot["Actual"] = self.actual_df.loc[mask]["Confirmed"].values
            
            # prepare plots
            # todo: make plots bigger
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
                # todo: what other variables should be visualized?
                _ = self.snl.history(target="Rt")
                _ = self.snl.history(target="rho").head()
                _ = self.snl.history(target="sigma").head()
                _ = self.snl.history(target="kappa").head()
                _ = self.snl.history(target="theta").head()

                
            # todo: how to stop describe from printing
            return self.infected_plot,self.confirmed_plot #,self.snl.describe()