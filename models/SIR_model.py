# -*- coding: utf-8 -*-
"""
Created on Thu May 12 19:50:49 2022

@author: User
"""

# imports
import pandas as pd
import covsirphy as cs
import matplotlib.pyplot as plt
import numpy as np
import sys
import warnings

# idea: implement __iter__, __str__
colors = ["r","g","b","m","black","brown"]


class SIR_model():

    def __init__(self, country, start_date, end_date):
        """
        Call function for setup.
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
        self.create=False

        self.color_list = colors
        self.colors = {"Main":"y","Actual":"c"}

    def load_data(self):
        """
        Loads up-to-date data from kaggle and stores it as object variable.

        """
        # data will be downloaded if existing files are older than 72h
        print("Checking if data is up-to-date... ")
        kaggle_path = sys.path[0] + '../kaggle/input'
        data_loader = cs.DataLoader(directory=kaggle_path, update_interval=72)

        # The number of cases and population values in jhu format
        print("Loading the data...")
        self.jhu_data = data_loader.jhu()
        # cases and deaths whole dataset
        self.total_df = self.jhu_data.total()
        # calculate start conditions given start date and cleaned data
        # days simulated
        self.total_days = (self.end_date - self.start_date).days

    def prep_data(self):
        """
        Extract clean data and subset for defined country.

        """
        # get rid of rows with missing data
        whole_df = self.jhu_data.cleaned()
        whole_country_df = whole_df.loc[(whole_df["Country"] == self.country)].copy()
        # change time format
        whole_country_df['Date'] = pd.to_datetime(whole_country_df['Date'])
        # store whole cleaned data
        self.whole_country_df = whole_country_df

        # drop provinces
        province_mask = self.whole_country_df["Province"] == "-"
        self.whole_country_df = self.whole_country_df.loc[province_mask]

        # define a time range
        mask = (self.whole_country_df['Date'] >= self.start_date) & (self.whole_country_df['Date'] <= self.end_date)
        df_timerange = self.whole_country_df.loc[mask]

        df_timerange.index = pd.RangeIndex(len(df_timerange))
        # store subsetted data
        self.df_timerange = df_timerange

        self.actual_df, _ = self.jhu_data.records(country="Switzerland")
        self.actual_df = self.actual_df.set_index("Date")

        # plotting the dataframe over the whole timeframe
        # self.actual_df.plot()

    def start_cond(self):
        """
        Extract starting conditions from data frame.

        """
        self.start_fatal_recovered = self.df_timerange["Fatal"][0] + self.df_timerange["Recovered"][0]
        self.start_infected = self.df_timerange["Infected"][0]
        self.start_pop = self.df_timerange["Population"][0]
        # as population size is constant, fatal cases have to be subtracted as well
        self.start_susceptible = self.start_pop - self.start_infected - self.start_fatal_recovered
        self.start_dict = {'Fatal or Recovered': self.start_fatal_recovered, 'Infected': self.start_infected,
                           'Susceptible': self.start_susceptible}


    def create_sir(self,params={"rho":None,"sigma":None},plot=False):

        """
        Creates SIR model with defined parameters.

        """
        self.create = True

        params = self.check_params(params)

        # set model parameters
        self.model = cs.SIR
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

    def check_params(self, params):
        """
        Check for correctness of input params.
        """
        
        if params["rho"] is None and params["sigma"] is None:
            print("No values for rho and sigma given. Estimation of both parameters, this may take a while...")
            params["rho"], params["sigma"] = self.parameter_estimation()
            print(f"Estimations: Rho = {params['rho']}, Sigma = {params['sigma']}\n\n")
            return params

        elif params["rho"] is None:
            print("No value for rho given. Estimation of rho, this may take a while...")
            params["rho"], _ = self.parameter_estimation()
            print(f"Estimation: Rho = {params['rho']}\n\n")
            return params

        elif params["sigma"] is None:
            print("No value for sigma given. Estimation of sigma, this may take a while...")
            _, params["sigma"] = self.parameter_estimation()
            print(f"Estimation: Sigma = {params['sigma']}\n\n")
            return params

        else:
            return params


    def create_main(self):
        """
        Computes values of SIR model without any changes of parameters for the time range given
        """
        
        if not self.create:
            raise Warning("create_sir has to be called before create_main")

        self.main = True
        self.scenario_names = []

        # use create_main as reset function for object
        self.snl = cs.Scenario(tau=1440, **self.area)
        self.snl.register(self.example_data)


        # if show_figure=True, this plots the main scenario
        record_df = self.snl.records(show_figure=False)

        # Set 0th phase from e.g. 01Sep2020 to 01Dec2020 with preset parameter values
        self.snl.clear(include_past=True)


        # past phase
        self.snl.add(end_date=self.end_date, model=self.model, **self.model.EXAMPLE["param_dict"])
        # store date of last scenario added to main
        self.last_phase_added=self.end_date
        
    
    def create_scenario(self, name, scenario_end_list,rho_constant_list=None,sigma_constant_list=None,plot=False):
        """
        Creating scenarios: gives ability to change parameters at specific points in times, which faciliates the modeling of real world measures.

        """


        # create main has to be called before create scenario
        if not self.main:
            raise Warning("create main has to be called before create scenario")
        if name not in self.scenario_names:
            self.scenario_names.append(name)
            self.colors[name] = self.color_list[len(self.scenario_names)-1]
        else:
            raise Warning("Scenario already exists")
        
        # Add new scenario
        self.snl.clear(name=name)
        

        for i,(phase_date, rho_constant, sigma_constant) in enumerate(zip(scenario_end_list,rho_constant_list,sigma_constant_list)):
            phase_date = pd.to_datetime(phase_date)
            # only add phase if it goes beyond last added phase
            if phase_date > self.last_phase_added:
                self.snl.add(end_date=phase_date, name="Main")
                self.last_phase_added = phase_date


            # always update rho and sigma. If no change is desired than same constant has to be given as input
            rho_new = self.snl.get("rho", phase="0th") * rho_constant

            sigma_new = self.snl.get("sigma", phase="0th") * sigma_constant

            # Add th i-th phase with the newly calculated params
            self.snl.add(end_date=phase_date, name=name, rho=rho_new,sigma=sigma_new)

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
        size_factor = 1.5  # sets the ylim value in relation to the maximum value of the actual data
        fig, (ax_inf, ax_conf) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(12, 12))
        ax_inf.set_title("Active Infected")
        ax_conf.set_title("Total confirmed cases")
        # plot actual and main infected
        ax_inf.plot(self.infected_plot.index,
                    self.infected_plot["Actual"],
                    label="Actual",
                    color=self.colors["Actual"],
                    alpha=0.5)
        ax_inf.plot(self.infected_plot.index,
                    self.infected_plot["Main"],
                    label="Main",
                    color=self.colors["Main"],
                    alpha=0.5)
        ax_inf.set_ylim(top=max(self.infected_plot["Actual"] * size_factor), bottom=0)

        # plot actual and main confirmed
        ax_conf.plot(self.confirmed_plot.index,
                     self.confirmed_plot["Actual"],
                     label="Actual",
                     color=self.colors["Actual"],
                     alpha=0.5)
        ax_conf.plot(self.confirmed_plot.index,
                     self.confirmed_plot["Main"],
                     label="Main",
                     color=self.colors["Main"],
                     alpha=0.5)
        ax_conf.set_ylim(top=max(self.confirmed_plot["Actual"] * size_factor), bottom=0)

        for name in self.scenario_names:
            ax_inf.plot(self.infected_plot.index,
                        self.infected_plot[name],
                        label=name,
                        color=self.colors[name],
                        alpha=0.5)
            ax_conf.plot(self.confirmed_plot.index,
                         self.confirmed_plot[name],
                         label=name,
                         color=self.colors[name],
                         alpha=0.5)
        plt.legend()

        if plot:
            rt_history = self.snl.history(target="Rt", show_figure=False)
            rho_history = self.snl.history(target="rho", show_figure=False)
            sigma_history = self.snl.history(target="sigma", show_figure=False)

            fig, (ax_rt, ax_rho, ax_sigma) = plt.subplots(nrows=3, ncols=1, sharex=True)

            temp_list = ["Main"]
            temp_list.extend(self.scenario_names)
            for name in temp_list:
                ax_rt.plot(rt_history.index, rt_history[name], label=name, color=self.colors[name], alpha=0.5)
                ax_rho.plot(rho_history.index, rho_history[name], label=name, color=self.colors[name], alpha=0.5)
                ax_sigma.plot(sigma_history.index, sigma_history[name], label=name, color=self.colors[name], alpha=0.5)
            plt.legend()
      
     
    def simulate_scenario(self, name):
        """
        SIR model representation of scenario specified with name
        """

        if name not in self.scenario_names:
            raise Warning("Scenario does not exist")
        df = self.snl.simulate(name = name)
        return df


    def get_res_df(self):
        """
        Getter-function for resulting datafram

        """
        try:
            return self.res_df

        except:
            raise Warning("No res_df. Run create_sir first")

    def parameter_estimation(self):
        """
        Analyses the phases of the wave of infections and estimates their respective parameters.

        Returns
        -------
        Estimation for rho and sigma
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

        estimation_df = srt.summary()

        # the mean of the first 3 phases will be used for a first estimation

        rho_est = round(estimation_df["rho"].iloc[[0, 1, 2]].mean(axis=0), 4)
        sigma_est = round(estimation_df["sigma"].iloc[[0, 1, 2]].mean(axis=0), 4)

        return rho_est, sigma_est



def one_scenario():
    """
    one possible example for the use our class
    """
    start_date = '2020-09-01'
    end_date = '2021-03-01'
    country = "Switzerland"
    a = SIR_model(country, start_date, end_date)
    a.create_sir(rho=0.1, sigma=0.08)
    a.create_main()
    a.create_scenario(name="Lockdown", scenario_end_list=["31Mar2021", "20Apr2021", "30Apr2021"],
                      rho_constant_list=[0.5, 2, 0.5], sigma_constant_list=[2, 0.5, 1], plot=True)


