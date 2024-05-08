"""main_old.py: Main file for waste heat potential calculator (ETH Zurich, MEST programme, case studies group 1, 2023/24)"""
__author__ = "Hanne Goericke and Florian Schubert"
__copyright__ = "Copyright 2024"
__credits__ = [""]
__license__ = ""
__version__ = "1.0.0"
__email__ = ""
__status__ = ""

# imports
import mapping
import pandas as pd
from datetime import datetime

# program options
sensitivity_analysis = True # True for sensitivity analysis and False for program execution with default parameter values
sensitivity_parameter = "spatial_radius" # parameter for sensitivity analysis (select from below), if above is True
# T_source_min, T_source_max, T_sink_min, T_sink_max, T_delta, P_source_min_single, P_sink_min, COP, spatial_radius

# NO NEED to change code below
# adapt input_heat_pump.xlsx for parameter values of heatpump and spatial combinations (including sensitivity analysis)



# paths
path_sources_unconventional = "./input_unconventional_heat_sources.xlsx"
path_source_datacenter = "./input_data_center.xlsx"
path_sinks_dh = "./input_district_heating.xlsx"
path_heatpump = "./input_heat_pump.xlsx"
path_sources_industries = "./input_industries.xlsx"
path_full_load_hours = "./input_full_load_hours.xlsx"
#path_europe_map = "./NUTS_RG_20M_2021_3035/NUTS_RG_20M_2021_3035.shx"



# terminal, file path and plot options / values
number_data_summary = 30
number_plot_combinations = 200

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
pd.set_option('display.max_rows', 200)

date_and_time = datetime.today().strftime('%Y-%m-%d_%H-%M-%S')



# read input files and create initial dataframes
sources, sinks, heatpump = mapping.read_input(path_sources_unconventional, path_source_datacenter,
                                              path_sources_industries, path_sinks_dh, path_heatpump,
                                              path_full_load_hours)
# TODO TESTING
#sources, sinks, heatpump = mapping.read_input_industries_only(path_sources_industries, path_sinks_dh, path_heatpump,
#                                                              path_full_load_hours)
number_sources_initial = len(sources)
number_sinks_initial = len(sinks)


# switch program execution between sensitivity analysis and default parameter values (without variation)

# sensitivity analysis:
if sensitivity_analysis:

    # prepare all parameter values and variation names for sensitivity analysis (and assign physical unit for output)

    sensitivity_parameter_values = heatpump.at[sensitivity_parameter, 'values_sensitivity']
    sensitivity_parameter_values.insert(0, heatpump.at[sensitivity_parameter, 'value'])
    sensitivity_parameter_values = list(set(sensitivity_parameter_values))
    sensitivity_parameter_values.sort()

    sensitivity_variable_name = sensitivity_parameter.replace('_','-')
    sensitivity_variable_unit = ""
    if (sensitivity_parameter.startswith("T")):
        sensitivity_variable_unit = "C"
    elif (sensitivity_parameter.startswith("P")):
        sensitivity_variable_unit = "MW"
    elif (sensitivity_parameter == "spatial_radius"):
        sensitivity_variable_unit = "m"
    scenario_names_plot = []
    scenario_names_path = []
    for value in sensitivity_parameter_values:
        scenario_names_plot.append(sensitivity_parameter + " = " + str(value) + " " + sensitivity_variable_unit)
        scenario_names_path.append(sensitivity_variable_name + "_" + str(value).replace('.','-')
                                   + sensitivity_variable_unit)



    print("[PROGRAM] Run sensitivity analysis for " + sensitivity_parameter + " = "
          + str(sensitivity_parameter_values) + " " + sensitivity_variable_unit)


    variation_target = "" # indicates whether variation of parameter values already affected sources, sinks or hasn't
                            # been executed so far



    # initial source filtering

    if (sensitivity_parameter == "T_source_min"):
        # sensitivity analysis variates T_source_min
        variation_sources = []
        for value in sensitivity_parameter_values:
            variation_sources.append(mapping.filter_single_temperature("sources", sources,
                                                                       value, heatpump.at['T_source_max', 'value']))
        sources = variation_sources
        variation_target = "source"
    elif (sensitivity_parameter == "T_source_max"):
        # sensitivity analysis variates T_source_max
        variation_sources = []
        for value in sensitivity_parameter_values:
            variation_sources.append(mapping.filter_single_temperature("sources", sources,
                                                                       heatpump.at['T_source_min', 'value'], value))
        sources = variation_sources
        variation_target = "source"
    else:
        # sensitivity analysis neither varies T_source_min, nor T_source_max
        sources = mapping.filter_single_temperature("sources", sources,
                                                    heatpump.at['T_source_min', 'value'],
                                                    heatpump.at['T_source_max', 'value'])

    if (variation_target == "source"):
        # iterate through sensitivity analysis variations if T_source_min or T_source_max is sensitivity parameter
        for i in range(len(sensitivity_parameter_values)):
            sources[i] = mapping.filter_single_power("sources", sources[i],
                                                     heatpump.at['P_source_min_single', 'value'])
    else:
        if (sensitivity_parameter == "P_source_min_single"):
            # sensitivity analysis variates P_source_min_single
            variation_sources = []
            for value in sensitivity_parameter_values:
                variation_sources.append(mapping.filter_single_power("sources", sources, value))
            sources = variation_sources
            variation_target = "source"
        else:
            # sensitivity analysis neither varies T_source_min, T_source_max, nor P_source_min_single
            sources = mapping.filter_single_power("sources", sources,
                                                  heatpump.at['P_source_min_single', 'value'])


    # initial sink filtering

    if (sensitivity_parameter == "T_sink_min"):
        # sensitivity analysis variates T_sink_min
        variation_sinks = []
        for value in sensitivity_parameter_values:
            variation_sinks.append(mapping.filter_single_temperature("sinks", sinks, value,
                                                                     heatpump.at['T_sink_max', 'value']))
        sinks = variation_sinks
        variation_target = "sink"
    elif (sensitivity_parameter == "T_sink_max"):
        # sensitivity analysis variates T_sink_max
        variation_sinks = []
        for value in sensitivity_parameter_values:
            variation_sinks.append(mapping.filter_single_temperature("sinks", sinks,
                                                                     heatpump.at['T_sink_min', 'value'], value))
        sinks = variation_sinks
        variation_target = "sink"
    else:
        # sensitivity analysis neither varies T_sink_min, nor T_sink_max
        sinks = mapping.filter_single_temperature("sinks", sinks,
                                                    heatpump.at['T_sink_min', 'value'],
                                                    heatpump.at['T_sink_max', 'value'])

    if (variation_target == "sink"):
        # iterate through sensitivity analysis variations if T_sink_min or T_sink_max is sensitivity parameter
        for i in range(len(sensitivity_parameter_values)):
            sinks[i] = mapping.filter_single_power("sinks", sinks[i],
                                                   heatpump.at['P_sink_min', 'value'])
    else:
        if (sensitivity_parameter == "P_sink_min"):
            # sensitivity analysis variates P_sink_min
            variation_sinks = []
            for value in sensitivity_parameter_values:
                variation_sinks.append(mapping.filter_single_power("sinks", sinks, value))
            sinks = variation_sinks
            variation_target = "sink"
        else:
            # sensitivity analysis neither varies T_sink_min, T_sink_max, nor P_sink_min
            sinks = mapping.filter_single_power("sinks", sinks,
                                                heatpump.at['P_sink_min', 'value'])



    # spatial combinations and filtering

    if (variation_target == "source"):
        # iterate through sensitivity analysis variations if T_source_min, T_source_max or P_source_min_single is sensitivity parameter
        variation_sinks = []
        variation_combinations = []
        plot_maps = []
        for i in range(len(sensitivity_parameter_values)):
            variation_sinks.append(sinks)
            combinations = mapping.combine_sinks_with_sources(sources[i], sinks, heatpump.at['spatial_radius', 'value'])
            combinations = mapping.filter_temperature_matching(sources[i], sinks, combinations, heatpump.at['T_delta', 'value'])
            combinations = mapping.filter_power_matching(sources[i], sinks, combinations,
                                                         heatpump.at['P_sink_min', 'value'], heatpump.at['COP', 'value'])
            combinations = mapping.score(combinations)
            variation_combinations.append(combinations)

            plot_maps.append(mapping.render_map(sources[i], sinks, combinations, number_plot_combinations))

        sinks = variation_sinks
        combinations = variation_combinations

    elif (variation_target == "sink"):
        # iterate through sensitivity analysis variations if T_sink_min, T_sink_max or P_sink_min is sensitivity parameter
        variation_sources = []
        variation_combinations = []
        plot_maps = []
        for i in range(len(sensitivity_parameter_values)):
            variation_sources.append(sources)
            combinations = mapping.combine_sinks_with_sources(sources, sinks[i], heatpump.at['spatial_radius', 'value'])
            combinations = mapping.filter_temperature_matching(sources, sinks[i], combinations, heatpump.at['T_delta', 'value'])
            if (sensitivity_parameter == "P_sink_min"):
                combinations = mapping.filter_power_matching(sources, sinks[i], combinations,
                                                             sensitivity_parameter_values[i], heatpump.at['COP', 'value'])
            else:
                combinations = mapping.filter_power_matching(sources, sinks[i], combinations,
                                                             heatpump.at['P_sink_min', 'value'],
                                                             heatpump.at['COP', 'value'])
            combinations = mapping.score(combinations)
            variation_combinations.append(combinations)

            plot_maps.append(mapping.render_map(sources, sinks[i], combinations, number_plot_combinations))

        sources = variation_sources
        combinations = variation_combinations

    else:
        variation_sources = []
        variation_sinks = []
        variation_combinations = []
        plot_maps = []

        if (sensitivity_parameter == "spatial_radius"):
            # sensitivity analysis variates spatial_radius
            for value in sensitivity_parameter_values:
                variation_sources.append(sources)
                variation_sinks.append(sinks)
                combinations = mapping.combine_sinks_with_sources(sources, sinks, value)
                combinations = mapping.filter_temperature_matching(sources, sinks, combinations,
                                                                   heatpump.at['T_delta', 'value'])
                combinations = mapping.filter_power_matching(sources, sinks, combinations,
                                                             heatpump.at['P_sink_min', 'value'],
                                                             heatpump.at['COP', 'value'])
                combinations = mapping.score(combinations)
                variation_combinations.append(combinations)

                plot_maps.append(mapping.render_map(sources, sinks, combinations, number_plot_combinations))


        elif (sensitivity_parameter == "T_delta"):
            combinations = mapping.combine_sinks_with_sources(sources, sinks, heatpump.at['spatial_radius', 'value'])
            combinations_tmp = combinations
            # sensitivity analysis variates T_delta
            # (while spatial_radius filter is identical for all variations and its result stored in combinations_tmp)
            for value in sensitivity_parameter_values:
                combinations = combinations_tmp
                variation_sources.append(sources)
                variation_sinks.append(sinks)
                combinations = mapping.filter_temperature_matching(sources, sinks, combinations, value)
                combinations = mapping.filter_power_matching(sources, sinks, combinations,
                                                             heatpump.at['P_sink_min', 'value'],
                                                             heatpump.at['COP', 'value'])
                combinations = mapping.score(combinations)
                variation_combinations.append(combinations)

                plot_maps.append(mapping.render_map(sources, sinks, combinations, number_plot_combinations))


        elif (sensitivity_parameter == "COP"):
            combinations = mapping.combine_sinks_with_sources(sources, sinks, heatpump.at['spatial_radius', 'value'])
            combinations = mapping.filter_temperature_matching(sources, sinks, combinations,
                                                               heatpump.at['T_delta', 'value'])
            combinations_tmp = combinations
            # sensitivity analysis variates COP
            # (while spatial_radius and T_delta filters are identical for all variations and their results stored in combinations_tmp)
            for value in sensitivity_parameter_values:
                combinations = combinations_tmp
                variation_sources.append(sources)
                variation_sinks.append(sinks)
                combinations = mapping.filter_power_matching(sources, sinks, combinations,
                                                             heatpump.at['P_sink_min', 'value'], value)
                combinations = mapping.score(combinations)
                variation_combinations.append(combinations)

                plot_maps.append(mapping.render_map(sources, sinks, combinations, number_plot_combinations))


        sources = variation_sources
        sinks = variation_sinks
        combinations = variation_combinations

    # plot and write output to files

    plot_number = mapping.render_plots(scenario_names_plot, sources, sinks, combinations,
                                       number_sources_initial, number_sinks_initial)

    mapping.write_output("sensitivity-analysis_" + sensitivity_parameter.replace('_', '-'),
                         scenario_names_path, sources, sinks, combinations, plot_maps, plot_number, heatpump,
                         date_and_time, number_data_summary)




# program execution with default parameter values (i.e., no sensitivity analysis):
else:
    print("[PROGRAM] Run calculations with default heatpump parameters")

    # filter sources
    sources = mapping.filter_single_temperature("sources", sources, heatpump.at['T_source_min', 'value'],
                                                 heatpump.at['T_source_max', 'value'])
    sources = mapping.filter_single_power("sources", sources, heatpump.at['P_source_min_single', 'value'])

    # filter sinks
    sinks = mapping.filter_single_temperature("sinks", sinks, heatpump.at['T_sink_min', 'value'],
                                             heatpump.at['T_sink_max', 'value'])
    sinks = mapping.filter_single_power("sinks", sinks, heatpump.at['P_sink_min', 'value'])

    # find and filter spatial combinations
    combinations = mapping.combine_sinks_with_sources(sources, sinks, heatpump.at['spatial_radius', 'value'])
    combinations = mapping.filter_temperature_matching(sources, sinks, combinations, heatpump.at['T_delta', 'value'])
    combinations = mapping.filter_power_matching(sources, sinks, combinations,
                                                 heatpump.at['P_sink_min', 'value'], heatpump.at['COP', 'value'])
    combinations = mapping.score(combinations)

    # create plots and write output to files

    plot_map = mapping.render_map(sources, sinks, combinations, number_plot_combinations)

    plot_number = mapping.render_plots(["filtered"], [sources], [sinks], [combinations],
                                                     number_sources_initial, number_sinks_initial)

    mapping.write_output("default-parameters", ["default-parameters"], [sources], [sinks],
                         [combinations], [plot_map], plot_number, heatpump, date_and_time,
                         number_data_summary)

print("[PROGRAM] Program finished")
