"""mapping.py: Mapping functions for waste heat potential calculator (ETH Zurich, MEST programme, case studies group 1, 2023/24)"""
__author__ = "Hanne Goericke and Florian Schubert"
__copyright__ = "Copyright 2024"
__credits__ = [""]
__license__ = ""
__version__ = "1.0.0"
__email__ = ""
__status__ = ""

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import matplotlib.pyplot as plt
import os

color_blue = '#1d344c'  # blue
color_gray = '#bcbbbc'  # gray
color_pink = '#e3144c'  # pink


# TODO print parameter mutation for sensible output in functions

def read_input(path_sources_unconventional, path_source_datacenter, path_sources_industries, path_sinks_dh,
               path_heatpump, path_full_load_hours):
    print("[DEBUG] load input data")

    # sinks datatypes:
    #   id                  int with source id (starting at 1, identical to index)
    #   latitude            float with latitude
    #   longitude           float with longitude
    #   name                string with plant name
    #   sector              string with sector name
    #   T_min               float with minimum temperature [°C]
    #   T_max               float with maximum temperature [°C]
    #   P                   float with power [MW]
    #   energy              float with energy [TJ]
    #   heatpumps           int[] with heatpump ids
    #   datasource          string with source name
    #   datasource_id       string with source id

    # Sink file:
    sinks = pd.read_excel(path_sinks_dh)
    sinks = sinks[sinks["latitude"] != "not geolocalized"]
    sinks = sinks[sinks["longitude"] != "not geolocalized"]
    sinks["latitude"] = sinks["latitude"].astype(float)
    sinks["longitude"] = sinks["longitude"].astype(float)
    sinks["heatpumps"] = [list() for x in range(len(sinks.index))]
    sinks["P"] = sinks["energy"] * (1 / (3600 * 8760) * 10 ** 6)
    sinks["datasource"] = "nan"
    sinks["datasource_id"] = "nan"
    sinks.index += 1

    # full load hours:
    full_load_hours = pd.read_excel(path_full_load_hours)
    full_load_hours.set_index('categories', inplace=True)
    # print(full_load_hours)

    # sources datatypes:
    #   id                  int with source id (starting at 1, identical to index)
    #   latitude            float with latitude
    #   longitude           float with longitude
    #   name                string with plant name
    #   sector              string with sector name
    #   T_min               float with minimum temperature [°C]
    #   T_max               float with maximum temperature [°C]
    #   P                   float with power [MW]
    #   energy              float with energy [TJ]
    #   heatpumps           int[] with heatpump ids
    #   datasource          string with source name
    #   datasource_id       string with source id

    # energy in GWh
    sources = pd.read_excel(path_sources_unconventional)
    sources = sources[sources["latitude"] != "not geolocalized"]
    sources = sources[sources["longitude"] != "not geolocalized"]
    sources["latitude"] = sources["latitude"].astype(float)
    sources["longitude"] = sources["longitude"].astype(float)
    sources["heatpumps"] = [list() for x in range(len(sources.index))]
    for s, sourcesunc in sources.iterrows():
        sources.at[s, "P"] = sourcesunc["energy"] * (
                    1 / (full_load_hours.at[str(sourcesunc["sector"]), 'full_load_hours']) * 10 ** 3)
    sources.drop(columns='energy', inplace=True)

    # power in MW
    sources_datacenter = pd.read_excel(path_source_datacenter)
    sources_datacenter = sources_datacenter[sources_datacenter["latitude"] != "not geolocalized"]
    sources_datacenter = sources_datacenter[sources_datacenter["longitude"] != "not geolocalized"]
    sources_datacenter["latitude"] = sources_datacenter["latitude"].astype(float)
    sources_datacenter["longitude"] = sources_datacenter["longitude"].astype(float)
    sources_datacenter["heatpumps"] = [list() for x in range(len(sources_datacenter.index))]

    # energy in TJ
    sources_industries = pd.read_excel(path_sources_industries)
    sources_industries["latitude"] = sources_industries["latitude"].astype(float)
    sources_industries["longitude"] = sources_industries["longitude"].astype(float)
    for s, sourcesind in sources_industries.iterrows():
        sources_industries.at[s, "P"] = sourcesind["energy"] * (
                    1 / (3600 * full_load_hours.at[str(sourcesind["sector"]), 'full_load_hours']) * 10 ** 6)

    sources = pd.concat([sources, sources_datacenter, sources_industries])
    sources['id'] = range(1, len(sources) + 1)
    sources.reset_index(drop=True, inplace=True)
    sources.index += 1

    # heatpump datatypes:
    #   property            string with property name (see below)
    #   value               float with default value for respective property
    #   values_sensitivity  float[] with alternative values for sensitivity analysis for respective property
    # properties:
    #   T_source_min            minimum source temperature [°C]
    #   T_source_max            maximum source temperature [°C]
    #   T_sink_min              minimum sink temperature [°C]
    #   T_sink_max              maximum sink temperature [°C]
    #   T_delta                 temperature difference for heat transfer [°C]
    #   P_source_min_single     minimum power for single source [MW]
    #   P_sink_min              minimum sink side (output) power [MW]
    #   COP                     average COP
    #   spatial_radius          allowed spatial radius [m]

    heatpump = pd.read_excel(path_heatpump, index_col='property')
    heatpump['values_sensitivity'] = heatpump['values_sensitivity'].apply(
        lambda x: [float(i) for i in x.split(',')] if isinstance(x, str) else [x])
    heatpump['values_sensitivity'] = heatpump['values_sensitivity'].apply(
        lambda x: [] if len(x) == 1 and pd.isna(x[0]) else x)

    return [sources, sinks, heatpump]


def read_input_industries_only(path_sources_industries, path_sinks_dh,
                               path_heatpump, path_full_load_hours):
    print("[DEBUG] load input data (industries only for testing purposes)")

    # sinks datatypes:
    #   id                  int with source id (starting at 1, identical to index)
    #   latitude            float with latitude
    #   longitude           float with longitude
    #   name                string with plant name
    #   sector              string with sector name
    #   T_min               float with minimum temperature [°C]
    #   T_max               float with maximum temperature [°C]
    #   P                   float with power [MW]
    #   energy              float with energy [TJ]
    #   datasource          string with source name
    #   datasource_id       string with source id

    # Sink file:
    sinks = pd.read_excel(path_sinks_dh)
    sinks = sinks[sinks["latitude"] != "not geolocalized"]
    sinks = sinks[sinks["longitude"] != "not geolocalized"]
    sinks["latitude"] = sinks["latitude"].astype(float)
    sinks["longitude"] = sinks["longitude"].astype(float)
    sinks["P"] = sinks["energy"] * (1 / (3600 * 8760) * 10 ** 6)
    sinks["datasource"] = "nan"
    sinks["datasource_id"] = "nan"
    sinks.index += 1

    # full load hours:
    full_load_hours = pd.read_excel(path_full_load_hours)
    full_load_hours.set_index('categories', inplace=True)

    # sources datatypes:
    #   id                  int with source id (starting at 1, identical to index)
    #   latitude            float with latitude
    #   longitude           float with longitude
    #   name                string with plant name
    #   sector              string with sector name
    #   T_min               float with minimum temperature [°C]
    #   T_max               float with maximum temperature [°C]
    #   P                   float with power [MW]
    #   energy              float with energy [TJ]
    #   datasource          string with source name
    #   datasource_id       string with source id

    # energy in TJ
    sources_industries = pd.read_excel(path_sources_industries)
    sources_industries["latitude"] = sources_industries["latitude"].astype(float)
    sources_industries["longitude"] = sources_industries["longitude"].astype(float)
    for s, sourcesind in sources_industries.iterrows():
        sources_industries.at[s, "P"] = sourcesind["energy"] * (
                1 / (3600 * full_load_hours.at[str(sourcesind["sector"]), 'full_load_hours']) * 10 ** 6)

    sources = sources_industries
    sources['id'] = range(1, len(sources) + 1)
    sources.reset_index(drop=True, inplace=True)
    sources.index += 1

    # heatpump datatypes:
    #   property            string with property name (see below)
    #   value               float with default value for respective property
    #   values_sensitivity  float[] with alternative values for sensitivity analysis for respective property
    # properties:
    #   T_source_min            minimum source temperature [°C]
    #   T_source_max            maximum source temperature [°C]
    #   T_sink_min              minimum sink temperature [°C]
    #   T_sink_max              maximum sink temperature [°C]
    #   T_delta                 temperature difference for heat transfer [°C]
    #   P_source_min_single     minimum power for single source [MW]
    #   P_sink_min              minimum sink side (output) power [MW]
    #   COP                     average COP
    #   spatial_radius          allowed spatial radius [m]

    heatpump = pd.read_excel(path_heatpump, index_col='property')
    heatpump['values_sensitivity'] = (heatpump['values_sensitivity'].apply(
        lambda x: [float(i) for i in x.split(',')] if isinstance(x, str) else [x]))
    heatpump['values_sensitivity'] = heatpump['values_sensitivity'].apply(
        lambda x: [] if len(x) == 1 and pd.isna(x[0]) else x)

    return [sources, sinks, heatpump]


def filter_single_temperature(dataframe_name, dataframe, minimum, maximum):
    print("[DEBUG] filter " + dataframe_name + " for heatpump temperature range "
          + "(with T_" + dataframe_name + "_min=" + str(minimum) + "°C, "
          + "T_" + dataframe_name + "_max=" + str(maximum) + "°C)")

    for id, row in dataframe.iterrows():
        if (row['T_min'] < minimum) or (maximum < row['T_max']):
            dataframe = dataframe.drop(index=id)

    return dataframe


def filter_single_power(dataframe_name, dataframe, minimum):
    print("[DEBUG] filter " + dataframe_name + " for minimum power "
          + "(with P_" + dataframe_name + "_min=" + str(minimum) + "MW)")

    for id, row in dataframe.iterrows():
        if (row['P'] < minimum):
            dataframe = dataframe.drop(index=id)

    return dataframe


def combine_sinks_with_sources(sources, sinks, radius):
    print("[DEBUG] combine sinks with sources according to spatial distribution "
          + "(with r=" + str(radius) + "m)")

    # Spatial radius (HP) is given in [m] and we need [°] --> 1° = 70 000m
    radius_degrees = radius / 70000

    # Generate geodataframes
    # sinks_geo = gpd.GeoDataFrame(sinks, geometry=gpd.points_from_xy(sinks.longitude, sinks.latitude))
    sources_geo = gpd.GeoDataFrame(sources, geometry=gpd.points_from_xy(sources.longitude, sources.latitude))

    # combinations datatypes:
    #   id                      int with combination id (starting at 1, identical to index)
    #   id_sink                 int with corresponding sink id
    #   id_sources              int[] with corresponding source ids
    #   T_matching_coefficient  float with temperature matching coefficient (lambda)
    #   P_matching_coefficient  float with power matching coefficient (theta)
    #   power_capacity          float with transferable power capacity [W]
    #   score                   float with score based on temperature/power coefficient and capacity (the higher, the better)
    combinations_columns = ['id', 'id_sink', 'id_sources', 'T_matching_coefficient',
                            'P_matching_coefficient', 'power_capacity', 'score']
    combinations_indices = []
    combinations_data = []
    id = 0
    # Iterate through each sink
    for sink_id, sink_row in sinks.iterrows():
        sink_point = Point(sink_row['longitude'], sink_row['latitude'])

        # Create buffer around sink point based on spatial radius
        buffer_zone = sink_point.buffer(radius_degrees)

        # Find sources within the buffer zone
        sources_within_buffer = sources_geo[sources_geo.geometry.within(buffer_zone)]

        # Collect all source ids within the buffer zone
        source_ids = sources_within_buffer['id'].tolist()

        # Check if source_ids list is not empty
        if source_ids:
            id = id + 1
            combinations_indices.append(id)
            combinations_data.append([id, sink_row['id'], source_ids, 0, 0, 0, 0])

    combinations = pd.DataFrame(columns=combinations_columns, index=combinations_indices, data=combinations_data)

    return combinations


def filter_temperature_matching(sources, sinks, combinations, T_delta):
    print("[DEBUG] filter spatial distribution map by temperatures and score temperature matching "
          + "(with T_delta=" + str(T_delta) + "°C)")

    for combination_id, combination_row in combinations.iterrows():
        sink_id = combination_row['id_sink']

        T_sink_min = sinks.at[sink_id, 'T_min']
        T_sink_max = sinks.at[sink_id, 'T_max']

        gamma_values = []
        source_power_total = 0

        filtered_sources = []
        for source_id in combination_row['id_sources']:
            T_source_min = sources.at[source_id, 'T_min']
            if (T_source_min + T_delta <= T_sink_max):
                filtered_sources.append(source_id)

        if (len(filtered_sources) > 0):
            combinations.at[combination_id, 'id_sources'] = tuple(filtered_sources)
            for source_id in combination_row['id_sources']:
                T_source_min = sources.at[source_id, 'T_min']
                T_source_max = sources.at[source_id, 'T_max']
                gamma_source = ((T_sink_min - T_delta - T_source_max) /
                                (T_sink_max - T_delta - T_source_min))
                gamma_values.append(gamma_source)
                source_power_total += sources.at[source_id, 'P']
            combinations.at[combination_id, 'T_matching_coefficient'] = sum(gamma_values) / source_power_total
        else:
            combinations = combinations.drop(index=combination_id)

    return combinations


def filter_power_matching(sources, sinks, combinations, P_sink_min, COP):
    print("[DEBUG] filter spatial distribution map by powers"
          + "(with P_sink_min=" + str(P_sink_min) + "MW, "
          + "COP=" + str(COP) + ")")

    for combination_id, combination_row in combinations.iterrows():
        capacity_sink = sinks.at[combination_row['id_sink'], 'P']
        capacity_sources = 0
        for source_id in combination_row['id_sources']:
            capacity_sources += sources.at[source_id, 'P']

        capacity_sources = capacity_sources * COP
        capacity = min(capacity_sink, capacity_sources)

        if (capacity < P_sink_min):
            combinations = combinations.drop(index=combination_id)
        else:
            combinations.at[combination_id, 'power_capacity'] = capacity
            combinations.at[combination_id, 'P_matching_coefficient'] = capacity / capacity_sink

    return combinations


def score(combinations):
    print("[DEBUG] score combinations")

    for combination_id, combination_row in combinations.iterrows():
        score = ((combination_row['T_matching_coefficient'] + combination_row['P_matching_coefficient']) *
                 combination_row['power_capacity'])
        combinations.at[combination_id, 'score'] = score

    spatial_combinations = combinations.sort_values(by=['score'], ascending=False)

    return spatial_combinations


def print_dataframes(sources, sinks, combinations, heatpump):
    print("----- SOURCES -----")
    print(sources)
    print("----- SINKS -----")
    print(sinks)
    print("----- SPATIAL COMBINATIONS -----")
    print(combinations)
    print("----- HEATPUMP -----")
    print(heatpump)


def render_map(sources, sinks, combinations, combination_number):

    sources = gpd.GeoDataFrame(sources, geometry=gpd.points_from_xy(sources.longitude, sources.latitude))
    sinks = gpd.GeoDataFrame(sinks, geometry=gpd.points_from_xy(sinks.longitude, sinks.latitude))

    path_to_europe = "NUTS_BN_20M_2021_4326/NUTS_BN_20M_2021_4326.shx"
    europe = gpd.read_file(path_to_europe)

    # Create figure with Europe map
    fig, ax = plt.subplots(figsize=(7, 7))

    # Plot Europe map
    europe.plot(color='grey', edgecolor='black', ax=ax, alpha=0.5, zorder=0)

    # Initialize dataframes
    df_sinks = pd.DataFrame(columns=sinks.columns)
    df_sources = pd.DataFrame(columns=sources.columns)

    # Iterate through spatial combinations
    for _, combination in combinations.head(combination_number).iterrows():
        # Find corresponding sink
        sink = sinks[sinks['id'] == combination['id_sink']].iloc[0]
        df_sinks = pd.concat([df_sinks, pd.DataFrame(sink).T])

        # Iterate through source IDs
        for source_id in combination['id_sources']:
            # Check if source ID exists in sources DataFrame
            if source_id in sources['id'].values:
                # Find corresponding source
                source = sources[sources['id'] == source_id].iloc[0]
                df_sources = pd.concat([df_sources, pd.DataFrame(source).T])

    # Convert 'P' column to float and adjust scale
    df_sinks['P'] = sinks['P'].astype(float)
    #df_sources['P'] = df_sources['P'].astype(float) *2


    # Plot sink and source scatter points
    ax.scatter(df_sinks['longitude'], df_sinks['latitude'], s=50*np.log(df_sinks['P']), label='Sinks', alpha=0.7, color=color_pink, zorder=1)
    #ax.scatter(df_sources['longitude'], df_sources['latitude'], s=df_sources['P'], label='Sources', alpha=1,
    #           color="red", zorder=1)

    # Set limits for x and y axes
    ax.set_xlim(-10, 30)
    ax.set_ylim(35, 70)

    # Hide the axes
    ax.set_axis_off()

    #plt.show()

    return fig




def render_plots(scenario_names, sources, sinks, combinations, initial_source_count, initial_sink_count):
    categories1 = ["input data (unfiltered)"]
    values_sources = [initial_source_count]
    values_sinks = [initial_sink_count]
    values_combinations = []
    for i in range(len(scenario_names)):
        categories1.append(scenario_names[i])
        values_sources.append(len(sources[i]))
        values_sinks.append(len(sinks[i]))
        values_combinations.append(len(combinations[i]))

    bar_width = 0.35

    ind1 = np.arange(len(categories1))
    fig1, ax1 = plt.subplots(figsize=(10, 2 + len(categories1)))
    ax1.barh(ind1 + bar_width / 1.9, values_sources, bar_width, color=color_blue, label='sources')
    ax1.barh(ind1 - bar_width / 1.9, values_sinks, bar_width, color=color_pink, label='sinks')
    ax1.set(yticks=ind1, yticklabels=categories1, xlabel='number of entries')
    ax1.legend()
    plt.tight_layout()

    if (len(scenario_names) > 1):
        categories2 = scenario_names
        ind2 = np.arange(len(categories2))
        fig2, ax2 = plt.subplots(figsize=(10, 2 + 0.5 * len(categories2)))
        ax2.barh(ind2, values_combinations, 2 * bar_width, color=color_gray, label='combinations')
        ax2.set(yticks=ind2, yticklabels=categories2, xlabel='number of entries')
        ax2.legend()
        plt.tight_layout()

        return [fig1, fig2]

    else:
        return [fig1]


def write_output(title, scenario_names, sources, sinks, combinations, plots_map, plots_sensitivity_analysis,
                 heatpump, date_and_time, combination_number):
    print("[DEBUG] write output data to files")

    path_directory = f"output_{date_and_time}_{title}"
    path_heatpump = f"{path_directory}/{title}_heatpump.xlsx"
    path_plot_sources_sinks = f"{path_directory}/{title}_plot_number-sinks-sources.pdf"
    path_plot_analysis_combinations = f"{path_directory}/{title}_plot_number-combinations.pdf"

    try:
        os.mkdir(path_directory)

        for i in range(len(scenario_names)):
            path_table = f"{path_directory}/{scenario_names[i]}_table.xlsx"
            path_summary = f"{path_directory}/{scenario_names[i]}_summary.txt"
            path_map = f"{path_directory}/{scenario_names[i]}_map.pdf"

            with pd.ExcelWriter(path_table, engine='xlsxwriter') as writer:
                combinations[i].to_excel(writer, sheet_name='Combinations', index=False)
                sinks[i].to_excel(writer, sheet_name='Sinks', index=False)
                sources[i].to_excel(writer, sheet_name='Sources', index=False)

            with open(path_summary, 'w') as f:
                counter = 0
                for combination_id, combination_row in combinations[i].iterrows():
                    if not (counter < combination_number and counter < len(combinations[i])):
                        break
                    counter += 1
                    f.write("----- ENTRY " + str(counter) + " -----\n\n")
                    f.write("combination:\n" + str(combination_row.to_frame().T) + "\n\n")
                    f.write("sink:\n" + str(sinks[i].loc[combination_row["id_sink"]].to_frame().T) + "\n\n")
                    f.write("sources:\n" + str(
                        sources[i].loc[sources[i]['id'].isin(combination_row["id_sources"])]) + "\n\n\n\n")

            plots_map[i].savefig(path_map)

        with pd.ExcelWriter(path_heatpump, engine='xlsxwriter') as writer:
            heatpump.to_excel(writer, sheet_name='Heatpump', index=True)

        plots_sensitivity_analysis[0].savefig(path_plot_sources_sinks)
        if (len(plots_sensitivity_analysis) > 1):
            plots_sensitivity_analysis[1].savefig(path_plot_analysis_combinations)

    except FileExistsError:
        print(f"Directory '{path_directory}' already exists. No data has been written.")
